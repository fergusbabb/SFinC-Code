import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt

lam = 1
gam = 1

def vectorPrime(state, N):
    global lam
    global gam
    
    x = state[0]
    y = state[1]

    xprime = (-3*x) + (lam * np.sqrt(1.5) * y**2)  + (1.5*x * (2*x**2 + gam*(1 - x**2 - y**2)))
    yprime =         (-lam * np.sqrt(1.5) * x * y) + (1.5*y * (2*x**2 + gam*(1 - x**2 - y**2)))

    return [xprime, yprime]

theta = np.linspace(0, np.pi, 51)
unitCircleX = np.cos(theta)
unitCircleY = np.sin(theta)

plt.plot(unitCircleX, unitCircleY, color="r")

N = np.linspace(0, 8, 201)

##x0 = -0.9
##y0 = 0.1
for x0 in np.arange(-0.99, 1, 0.1):
    for y0 in [0.01]:#np.arange(0, 1, 0.1):
        if x0**2 + y0**2 > 1:
            continue
        path = integrate.odeint(vectorPrime, [x0, y0], N)
        pathTranspose = path.transpose()

        pathx = pathTranspose[0]
        pathy = pathTranspose[1]

        plt.plot(pathx, pathy, color="b")

plt.axis("equal")
plt.show()
##plt.plot()
