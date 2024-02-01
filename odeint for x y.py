import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt

lam = 1
gam = 1

plt.annotate("$\lambda = $" + str(lam) + "\n$\gamma = $" + str(gam), xy = (-1, 1))

def vectorPrime(state, N):
    global lam
    global gam
    
    x = state[0]
    y = state[1]

    xprime = (-3*x) + (lam * np.sqrt(1.5) * y**2)  + (1.5*x * (2*x**2 + gam*(1 - x**2 - y**2)))
    yprime =         (-lam * np.sqrt(1.5) * x * y) + (1.5*y * (2*x**2 + gam*(1 - x**2 - y**2)))

    return [xprime, yprime]


if lam**2 > 3*gam:
    x = np.sqrt(1.5) * gam/lam
    y = np.sqrt(1.5 * (2 - gam) * gam/lam**2)

elif lam**2 < 6:
    x = lam/np.sqrt(6)
    y = np.sqrt(1 - x**2)

else:
    x = 1
    y = 0

omegaDE = x**2 + y**2


theta = np.linspace(0, np.pi, 51)
unitCircleX = np.cos(theta)
unitCircleY = np.sin(theta)

omegaDEcircleX = np.sqrt(omegaDE) * unitCircleX
omegaDEcircleY = np.sqrt(omegaDE) * unitCircleY

plt.plot(unitCircleX, unitCircleY, color="r")
plt.plot(omegaDEcircleX, omegaDEcircleY, "r--", linewidth=0.75)
plt.plot([0, x/np.sqrt(omegaDE)], [0, y/np.sqrt(omegaDE)], "r--", linewidth=0.75)

N = np.linspace(0, 8, 401)

for x0 in np.arange(-0.99, 1, 0.2):
    for y0 in np.arange(0.01, 1, 0.1):
        if x0**2 + y0**2 > 1:
            continue
        path = integrate.odeint(vectorPrime, [x0, y0], N)
        pathTranspose = path.transpose()

        pathx = pathTranspose[0]
        pathy = pathTranspose[1]

        plt.plot(pathx, pathy, color="k", linewidth=0.25)

plt.axis("equal")
plt.show()
