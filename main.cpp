//
// Created by a_d_e on 17.04.2021.
//
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

#define         HIGH_PRECISION    0.5
#define         LOW_PRECISION     1.0
#define         ZOOM_PRECISION    5.0
#define         ARROW_START       4
#define         MAX_SPLINE_STEP   0.01

#define Q(s)  (-(s))
#define EQN_NUMERATOR(dim) \
  (p0.dim*A_blend[0]+p1.dim*A_blend[1]+p2.dim*A_blend[2]+p3.dim*A_blend[3])

class Point {
public:
    double x, y;

    Point(double x, double y) : x(x), y(y) {}

    Point operator*(double a) {
        x *= a;
        y *= a;
        return *this;
    }

    Point operator+(Point point) {
        x += point.x;
        y += point.y;
        return *this;
    }

    Point operator/(double a) {
        x /= a;
        y /= a;
        return *this;
    }

    friend ostream &operator<<(ostream &os, const Point &point) {
        os << "x: " << point.x << "," << " y: " << point.y;
        return os;
    }
};

double f_blend(double numerator, double denominator) {
    double p = 2 * denominator * denominator;
    double u = numerator / denominator;
    double u2 = u * u;

    return (u * u2 * (10 - p + (2 * p - 15) * u + (6 - p) * u2));
}

double g_blend(double u, double q) {/* p equals 2 */
    return (u * (q + u * (2 * q + u * (8 - 12 * q + u * (14 * q - 11 + u * (4 - 5 * q))))));
}

double h_blend(double u, double q) {
    double u2 = u * u;
    return (u * (q + u * (2 * q + u2 * (-2 * q - u * q))));
}

void negative_s1_influence(double t, double s1, double *A0, double *A2) {
    *A0 = h_blend(-t, Q(s1));
    *A2 = g_blend(t, Q(s1));
}

void negative_s2_influence(double t, double s2, double *A1, double *A3) {
    *A1 = g_blend(1 - t, Q(s2));
    *A3 = h_blend(t - 1, Q(s2));
}

void positive_s1_influence(int k, double t, double s1, double *A0, double *A2) {
    double Tk;

    Tk = k + 1 + s1;
    *A0 = (t + k + 1 < Tk) ? f_blend(t + k + 1 - Tk, k - Tk) : 0.0;

    Tk = k + 1 - s1;
    *A2 = f_blend(t + k + 1 - Tk, k + 2 - Tk);
}

void positive_s2_influence(int k, double t, double s2, double *A1, double *A3) {
    double Tk;

    Tk = k + 2 + s2;
    *A1 = f_blend(t + k + 1 - Tk, k + 1 - Tk);

    Tk = k + 2 - s2;
    *A3 = (t + k + 1 > Tk) ? f_blend(t + k + 1 - Tk, k + 3 - Tk) : 0.0;
}

void add_point(Point point) {
    cout << point << endl;
}

void point_adding(double *A_blend, Point p0, Point p1, Point p2, Point p3) {
    double weights_sum;

    weights_sum = A_blend[0] + A_blend[1] + A_blend[2] + A_blend[3];
    add_point(Point(EQN_NUMERATOR(x) / (weights_sum), EQN_NUMERATOR(y) / (weights_sum)));
}

void point_computing(double *A_blend, Point p0, Point p1, Point p2, Point p3, double *x, double *y) {
    double weights_sum;

    weights_sum = A_blend[0] + A_blend[1] + A_blend[2] + A_blend[3];

    *x = EQN_NUMERATOR(x) / (weights_sum);
    *y = EQN_NUMERATOR(y) / (weights_sum);
}

float step_computing(int k, Point p0, Point p1, Point p2, Point p3, double s1, double s2, float precision) {
    double A_blend[4];
    double xstart, ystart, xend, yend, xmid, ymid, xlength, ylength;
    int start_to_end_dist, number_of_steps;
    float step, angle_cos, scal_prod, xv1, xv2, yv1, yv2, sides_length_prod;

    /* This function computes the step used to draw the segment (p1, p2)
       (xv1, yv1) : coordinates of the vector from middle to origin
       (xv2, yv2) : coordinates of the vector from middle to extremity */

    if ((s1 == 0) && (s2 == 0))
        return (1.0);              /* only one step in case of linear segment */

    /* compute coordinates of the origin */
    if (s1 > 0) {
        if (s2 < 0) {
            positive_s1_influence(k, 0.0, s1, &A_blend[0], &A_blend[2]);
            negative_s2_influence(0.0, s2, &A_blend[1], &A_blend[3]);
        } else {
            positive_s1_influence(k, 0.0, s1, &A_blend[0], &A_blend[2]);
            positive_s2_influence(k, 0.0, s2, &A_blend[1], &A_blend[3]);
        }
        point_computing(A_blend, p0, p1, p2, p3, &xstart, &ystart);
    } else {
        xstart = p1.x;
        ystart = p1.y;
    }

    /* compute coordinates  of the extremity */
    if (s2 > 0) {
        if (s1 < 0) {
            negative_s1_influence(1.0, s1, &A_blend[0], &A_blend[2]);
            positive_s2_influence(k, 1.0, s2, &A_blend[1], &A_blend[3]);
        } else {
            positive_s1_influence(k, 1.0, s1, &A_blend[0], &A_blend[2]);
            positive_s2_influence(k, 1.0, s2, &A_blend[1], &A_blend[3]);
        }
        point_computing(A_blend, p0, p1, p2, p3, &xend, &yend);
    } else {
        xend = p2.x;
        yend = p2.y;
    }

    /* compute coordinates  of the middle */
    if (s2 > 0) {
        if (s1 < 0) {
            negative_s1_influence(0.5, s1, &A_blend[0], &A_blend[2]);
            positive_s2_influence(k, 0.5, s2, &A_blend[1], &A_blend[3]);
        } else {
            positive_s1_influence(k, 0.5, s1, &A_blend[0], &A_blend[2]);
            positive_s2_influence(k, 0.5, s2, &A_blend[1], &A_blend[3]);
        }
    } else if (s1 < 0) {
        negative_s1_influence(0.5, s1, &A_blend[0], &A_blend[2]);
        negative_s2_influence(0.5, s2, &A_blend[1], &A_blend[3]);
    } else {
        positive_s1_influence(k, 0.5, s1, &A_blend[0], &A_blend[2]);
        negative_s2_influence(0.5, s2, &A_blend[1], &A_blend[3]);
    }
    point_computing(A_blend, p0, p1, p2, p3, &xmid, &ymid);

    xv1 = xstart - xmid;
    yv1 = ystart - ymid;
    xv2 = xend - xmid;
    yv2 = yend - ymid;

    scal_prod = xv1 * xv2 + yv1 * yv2;

    sides_length_prod = sqrt((xv1 * xv1 + yv1 * yv1) * (xv2 * xv2 + yv2 * yv2));

    /* compute cosinus of origin-middle-extremity angle, which approximates the
       curve of the spline segment */
    if (sides_length_prod == 0.0)
        angle_cos = 0.0;
    else
        angle_cos = scal_prod / sides_length_prod;

    xlength = xend - xstart;
    ylength = yend - ystart;

    start_to_end_dist = sqrt((double) xlength * (double) xlength + (double) ylength * (double) ylength);

    /* more steps if segment's origin and extremity are remote */
    number_of_steps = sqrt(start_to_end_dist) / 2;

    /* more steps if the curve is high */
    number_of_steps += (int) ((1 + angle_cos) * 10);

    if (number_of_steps == 0)
        step = 1;
    else
        step = precision / number_of_steps;

    if ((step > MAX_SPLINE_STEP) || (step == 0))
        step = MAX_SPLINE_STEP;
    return (step);
}

void spline_segment_computing(float step, int k, Point p0, Point p1, Point p2, Point p3, double s1, double s2) {
    double A_blend[4];

    if (s1 < 0) {
        if (s2 < 0) {
            for (double t = 0.0; t < 1; t += step) {
                negative_s1_influence(t, s1, &A_blend[0], &A_blend[2]);
                negative_s2_influence(t, s2, &A_blend[1], &A_blend[3]);

                point_adding(A_blend, p0, p1, p2, p3);
            }
        } else {
            for (double t = 0.0; t < 1; t += step) {
                negative_s1_influence(t, s1, &A_blend[0], &A_blend[2]);
                positive_s2_influence(k, t, s2, &A_blend[1], &A_blend[3]);

                point_adding(A_blend, p0, p1, p2, p3);
            }
        }
    } else if (s2 < 0) {
        for (double t = 0.0; t < 1; t += step) {
            positive_s1_influence(k, t, s1, &A_blend[0], &A_blend[2]);
            negative_s2_influence(t, s2, &A_blend[1], &A_blend[3]);

            point_adding(A_blend, p0, p1, p2, p3);
        }
    } else {
        for (double t = 0.0; t < 1; t += step) {
            positive_s1_influence(k, t, s1, &A_blend[0], &A_blend[2]);
            positive_s2_influence(k, t, s2, &A_blend[1], &A_blend[3]);

            point_adding(A_blend, p0, p1, p2, p3);
        }
    }
}


bool compute_open_spline(vector<Point> points, double *sfactors, float precision) {
    float step;
    for (int i = 0; i < points.size() - 2; i++) {
        step = step_computing(i, points[i], points[i + 1], points[i + 2], points[i + 3], sfactors[i + 1],
                              sfactors[i + 2], precision);
        spline_segment_computing(step, i, points[i], points[i + 1], points[i + 2], points[i + 3], sfactors[i + 1],
                                 sfactors[i + 2]);
    }
    return true;
}

int main() {
    vector<Point> points;
    vector<double> x = {1, 1, 3, 5, 7, 9, 9};
    vector<double> y = {4, 4, 6, 4, 6, 4, 4};
    double s[] = {0, 0, 1, 1, 1, 1, 0};
    for (int i = 0; i < x.size(); ++i) {
        Point point(x[i], y[i]);
        points.push_back(point);
//        cout << point << endl;
    }
    compute_open_spline(points, s, 0.001);
//    Xspline xspline(points, s);
//    ofstream file;
//    file.open("result.txt");
//    for (int i = 0; i < points.size() - 3; ++i) {
//        for (double t = 0; t <= 1; t += 0.01) {
//            file << setprecision(6) << fixed << xspline.C(t, i) << endl;
//        }
//    }
//    file.close();
    return 0;
}
