//
// Created by a_d_e on 17.04.2021.
//
#include <iostream>
#include <vector>

using namespace std;

class Point {
public:
    double x, y;

    Point(double x, double y) : x(x), y(y) {}

    Point operator*(double a) const {
        x * a;
        y * a;
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
        os << "x: " << point.x << " y: " << point.y;
        return os;
    }
};

class Xspline {
public:
    Xspline(const vector<Point> &points, double *s) : points(points), s(s) {
    }

    Point C(double t, int k) {
        return (points[k] * a0(t, s[k + 1]) + points[k + 1] * a1(t, s[k + 2]) + points[k + 2] * a2(t, s[k + 1]) +
                points[k + 3] * a3(t, s[k + 2])) / (a0(t, s[k + 1]) + a1(t, s[k + 2]) + a2(t, s[k + 1]) +
                                                    a3(t, s[k + 2]));
    }

private:
    double *s;
    vector<Point> points;

    double F(double u, double p) {
        return u * u * u * (10 - p + (2 * p - 15) * u + (6 - p) * u * u);
    }

    double f(double n, double d) {
        return F(n / d, 2 * d * d);
    }

    double g(double u, double q) {
        return u * (q + u * (2 * q + u * (8 - 12 * q + u * (14 * q - 11 + u * (4 - 5 * q)))));
    }

    double h(double u, double q) {
        return u * (q + u * (2 * q + u * u * (-2 * q - u * q)));
    }


    double a0(double t, double s) {
        if (s < 0) {
            return h(-t, -s);
        }
        if (s >= 0 and t < s) {
            return f(t - s, -1 - s);
        } else {
            return 0;
        }
    }

    double a1(double t, double s) {
        if (s < 0) {
            return g(1 - t, -s);
        } else {
            return f(t - 1 - s, -1 - s);
        }
    }

    double a2(double t, double s) {
        if (s < 0) {
            return g(t, -s);
        } else {
            return f(t + s, 1 + s);
        }
    }

    double a3(double t, double s) {
        if (s < 0) {
            return h(t - 1, -s);
        }
        if (s >= 0 and t > 1 - s) {
            return f(t - 1 + s, 1 + s);
        } else {
            return 0;
        }
    }

};

int main() {
    vector<Point> points;
    double x[] = {.1, .3, .5, .7, .9};
    double y[] = {.4, .6, .4, .6, .4};
    double s[] = {1, 1, 1, 1, 1, 1, 1};
    for (int i = 0; i < 7; ++i) {
        Point point(x[i], y[i]);
        points.push_back(point);
    }
    Xspline xspline(points, s);
    for (int i = 0; i < 3; ++i) {
        cout << xspline.C(0, i) << endl;
    }

    return 0;
}