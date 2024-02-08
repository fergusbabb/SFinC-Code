import numpy as np
from matplotlib import pyplot as plt

gam = 0.5

def acceleration(x, y):
    return -3*x**2 - 1.5*gam*(1 - x**2 - y**2) + 1

def acceleration_y_boundary(x):
    ysquared = (2/gam - 1) * x**2 - 2/(3*gam) + 1#, np.zeros(len(x)))
    ysquared = np.where(ysquared>0, ysquared, 0)
    return np.sqrt(ysquared)

x = np.linspace(-1, 1, 401)

theta = np.linspace(0, np.pi, 51)
unitCircleX = np.cos(theta)
unitCircleY = np.sin(theta)

xWhere = np.where(acceleration_y_boundary(x) < np.sqrt(1-x**2), True, False)
plt.plot(unitCircleX, unitCircleY, color="r")
plt.fill_between(x, acceleration_y_boundary(x), np.sqrt(1-x**2), where = xWhere)

plt.axis("equal")
plt.show()
