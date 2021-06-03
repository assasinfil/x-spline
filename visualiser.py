import matplotlib.pyplot as plt

while True:
    with open("cmake-build-debug/result.txt") as f:
        data = f.readlines()

    points_x = list()
    points_y = list()
    data.sort()
    for line in data:
        x, y = line.rstrip().split(', ')
        points_x.append(float(x[3:]))
        points_y.append(float(y[3:]))

    plt.plot(points_x, points_y)
    plt.show()
