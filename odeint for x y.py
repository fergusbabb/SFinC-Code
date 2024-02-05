import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt

lam = 3
gam = 1

fig, ax = plt.subplots(2, 2)

ax[0, 0].annotate("$\lambda = $" + str(lam) + "\n$\gamma = $" + str(gam), xy = (-1, 0.8))

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

ax[0, 0].plot(unitCircleX, unitCircleY, color="r")
ax[0, 0].plot(omegaDEcircleX, omegaDEcircleY, "r--", linewidth=0.75)
ax[0, 0].plot([0, x/np.sqrt(omegaDE)], [0, y/np.sqrt(omegaDE)], "r--", linewidth=0.75)

N = np.linspace(0, 8, 401)

for x0 in np.linspace(-0.99, 0.99, 10):
    for y0 in [0.01]:#np.arange(0.01, 1, 0.1):
        if x0**2 + y0**2 > 1:
            continue
        path = integrate.odeint(vectorPrime, [x0, y0], N)
        pathTranspose = path.transpose()

        pathx = pathTranspose[0]
        pathy = pathTranspose[1]

        ax[0, 0].plot(pathx, pathy, linewidth=0.25)

        ax[1, 0].plot(N, pathx)
        ax[1, 1].plot(N, pathy)

        ax[0, 1].plot(N, 1 - 3*pathx**2 + 1.5*gam*(1 - pathx**2 - pathy**2))

ax[0, 0].axis("equal")
plt.show()
