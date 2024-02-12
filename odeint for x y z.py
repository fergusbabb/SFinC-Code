import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt

lam = 1

fig, ax = plt.subplots(2, 2)

ax[0, 0].annotate("$\lambda = $" + str(lam), xy = (-1, 0.8))

def vectorPrime(state, N):
    global lam
    
    x = state[0]
    y = state[1]
    z = state[2]

    HdotSection = 1 + x**2 - y**2 - (1/3)*z**2

    xprime = (-3*x) + (lam * np.sqrt(1.5) * y**2)  + (1.5*x * HdotSection)
    yprime =         (-lam * np.sqrt(1.5) * x * y) + (1.5*y * HdotSection)
    zprime = (-2*z)                                + (1.5*z * HdotSection)

    return [xprime, yprime, zprime]

def accelerationExpression(x, y, z):
    HdotSection = 1 + x**2 - y**2 - (1/3)*z**2
    return -3/2 * HdotSection + 1    

theta = np.linspace(0, np.pi, 51)
unitCircleX = np.cos(theta)
unitCircleY = np.sin(theta)



ax[0, 0].plot(unitCircleX, unitCircleY, color="r")

N = np.linspace(0, 8, 401)
a = np.exp(N)


z0 = 0.98

y0 = 0.1
x0 = 0.1

print(x0**2 + y0**2 + z0**2)

path = integrate.odeint(vectorPrime, [x0, y0, z0], N)
pathTranspose = path.transpose()

pathx = pathTranspose[0]
pathy = pathTranspose[1]
pathz = pathTranspose[2]

ax[0, 0].plot(pathx, pathy, linewidth=0.75)

ax[1, 0].plot(N, pathx**2 + pathy**2)

ax[1, 1].plot(N, pathz**2,                           "r", label = "$\Omega_r = z^2$")
ax[1, 1].plot(N, 1 - pathx**2 - pathy**2 - pathz**2, "g", label = "$\Omega_m = 1 - x^2 - y^2 - z^2$")
ax[1, 1].plot(N, pathx**2 + pathy**2,                "b", label = "$\Omega_\phi = x^2 + y^2$")

ax[0, 1].plot(N, accelerationExpression(pathx, pathy, pathz))
##        ax[0, 1].plot(N, -3*pathx**2 - 1.5*gam*(1 - pathx**2 - pathy**2) + 1)

ax[0, 0].axis("equal")
ax[0, 0].set_xlabel("x")
ax[0, 0].set_ylabel("y")

ax[0, 1].set_xlabel("N")
ax[0, 1].set_ylabel("Acceleration")

ax[1, 0].set_xlabel("N")
ax[1, 0].set_ylabel("$x^2 + y^2 = \Omega_\phi$")

ax[1, 1].set_xlabel("N")
ax[1, 1].set_ylabel("Density Parameters")
ax[1, 1].legend()
plt.show()
