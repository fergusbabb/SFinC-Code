import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt

lam = 2 * np.sqrt(2/3)# - 0.2

fig, ax = plt.subplots(2, 2)

ax[0, 0].set_axis_off()
ax[0, 0] = fig.add_subplot(2, 2, 1, projection="3d")

##ax[0, 0].annotate("$\lambda = $" + str(lam), xy = (-1, 0.8))

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
    return -(3/2) * HdotSection + 1    

theta = np.linspace(0, np.pi, 51)
unitCircleX = np.cos(theta)
unitCircleY = np.sin(theta)

ax[0, 0].plot(unitCircleX, unitCircleY, [0 for x in unitCircleX], color="r")
ax[0, 0].plot(unitCircleX, [0 for x in unitCircleX], unitCircleY, color="r")
ax[0, 0].plot([0 for x in unitCircleX[:26]], unitCircleX[:26], unitCircleY[:26], color="r")

ax[0, 0].plot([-1, 1], [0, 0], [0, 0], "r")

N = np.linspace(0, 10, 401)
a = np.exp(N)

ax[0, 1].plot([0, N[-1]], [0, 0], "k--")

z0 = 0.9

y0 = 0.1
x0 = 0.1

print(x0**2 + y0**2 + z0**2)

path = integrate.odeint(vectorPrime, [x0, y0, z0], N)
pathTranspose = path.transpose()

pathx = pathTranspose[0]
pathy = pathTranspose[1]
pathz = pathTranspose[2]

ax[0, 0].plot(pathx, pathy, pathz, linewidth=0.75)

ax[1, 0].plot(N, pathx**2 + pathy**2)

ax[1, 1].plot(N, pathz**2,                           "r", label = "$\Omega_r = z^2$")
ax[1, 1].plot(N, 1 - pathx**2 - pathy**2 - pathz**2, "g", label = "$\Omega_m = 1 - x^2 - y^2 - z^2$")
ax[1, 1].plot(N, pathx**2 + pathy**2,                "b", label = "$\Omega_\phi = x^2 + y^2$")

ax[0, 1].plot(N, accelerationExpression(pathx, pathy, pathz))


##xFixed = 2 * np.sqrt(6) / lam
##yFixed = 2 / (np.sqrt(3) * lam)
##
##ax[0, 0].scatter(xFixed, yFixed, 0)

##ax[0, 0].axis("equal")
ax[0, 0].set_xlabel("x")
ax[0, 0].set_ylabel("y")
ax[0, 0].set_zlabel("z")
ax[0, 0].set_box_aspect([2, 1, 1])
ax[0, 0].set_xticks([-1, -0.5, 0, 0.5, 1])
ax[0, 0].set_yticks([0, 0.5, 1])
ax[0, 0].set_zticks([0, 0.5, 1])

ax[0, 1].set_xlabel("N")
ax[0, 1].set_ylabel("Acceleration")

ax[1, 0].set_xlabel("N")
ax[1, 0].set_ylabel("$x^2 + y^2 = \Omega_\phi$")

ax[1, 1].set_xlabel("N")
ax[1, 1].set_ylabel("Density Parameters")
ax[1, 1].legend()
plt.show()
