import matplotlib.pyplot as plt

with open("result.txt") as f:
    data = f.readlines()

points_x = list()
points_y = list()

for line in data:
    x, y = line.rstrip().split(', ')
    points_x.append(float(x[3:]))
    points_y.append(float(y[3:]))

plt.plot(points_x, points_y)
plt.show()
