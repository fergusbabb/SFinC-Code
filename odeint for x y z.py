import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt

lam = 2# np.sqrt(12/3)+0.1#2 * np.sqrt(2/3)# - 0.2

fig, ax = plt.subplots(2, 2)

ax[0, 0].set_axis_off()
ax[0, 0] = fig.add_subplot(2, 2, 1, projection="3d")


def vectorPrime(state, N):
    global lam
    
    x = state[0]
    y = state[1]
    z = state[2]

    HdotSection = 1 + x**2 - y**2 + (1/3)*z**2

    xprime = (-3*x) + (lam * np.sqrt(1.5) * y**2)  + (1.5*x * HdotSection)
    yprime =         (-lam * np.sqrt(1.5) * x * y) + (1.5*y * HdotSection)
    zprime = (-2*z)                                + (1.5*z * HdotSection)

    return [xprime, yprime, zprime]

def accelerationExpression(x, y, z):
    HdotSection = 1 + x**2 - y**2 + (1/3)*z**2
    return (1/2) * (-3*x**2 + 3*y**2 - z**2 - 1)

def gamma_phi(x, y):
    return (2*x**2) / (x**2 + y**2)

theta = np.linspace(0, np.pi, 51)
unitCircleX = np.cos(theta)
unitCircleY = np.sin(theta)

ax[0, 0].plot(unitCircleX, unitCircleY, [0 for x in unitCircleX], color="r")
ax[0, 0].plot(unitCircleX, [0 for x in unitCircleX], unitCircleY, color="r")
ax[0, 0].plot([0 for x in unitCircleX[:26]], unitCircleX[:26], unitCircleY[:26], color="r")

ax[0, 0].plot([-1, 1], [0, 0], [0, 0], "r")

N = np.linspace(0, 30, 401)
reverseN = -N
a = np.exp(N)

ax[0, 1].plot([0, N[-1]], [0, 0], "k--")

ax[1, 0].plot([0, N[-1]], [4/3, 4/3], "k--", linewidth = 0.5)
ax[1, 0].plot([0, N[-1]], [1, 1], "k--", linewidth = 0.5)


## Initial conditions

if False:#lam > 2:
    y0 = 2 / (lam * np.sqrt(3))#0.2
    x0 = 2 * np.sqrt(2/3) / lam #np.sqrt(2) * y0
    z0 = np.sqrt(1 - (4/lam**2)) - 0.00000005
else:
    x0 = 0#-0.01#0.05
    y0 = 0.01
    z0 = 0.9

print(x0**2 + y0**2 + z0**2)

fixedPoints = np.array([
    [0, 0, 0],
    [0, 0, 1],
    [1, 0, 0],
    [-1, 0, 0]
    ])

if lam**2 < 6:
    fixedPoints = np.append(fixedPoints, [[lam/np.sqrt(6), np.sqrt(1 - (lam**2)/6), 0]], axis=0)

if lam**2 > 3:
    fixedPoints = np.append(fixedPoints, [[np.sqrt(3/2)/lam, np.sqrt(3/2)/lam, 0]], axis=0)

if lam**2 > 8/3:
    fixedPoints = np.append(fixedPoints, [[2 * np.sqrt(2/3) / lam, 2 / (lam * np.sqrt(3)), np.sqrt(1 - (4/lam**2))]], axis=0)

path = integrate.odeint(vectorPrime, [x0, y0, z0], N)
pathTranspose = path.transpose()

pathx = pathTranspose[0]
pathy = pathTranspose[1]
pathz = pathTranspose[2]

ax[0, 0].plot(pathx, pathy, pathz, linewidth=0.75)

ax[1, 0].plot(N, gamma_phi(pathx, pathy))

ax[1, 1].plot(N, pathz**2,                           "r", label = "$\Omega_r = z^2$")
ax[1, 1].plot(N, 1 - pathx**2 - pathy**2 - pathz**2, "g", label = "$\Omega_m = 1 - x^2 - y^2 - z^2$")
ax[1, 1].plot(N, pathx**2 + pathy**2,                "b", label = "$\Omega_\phi = x^2 + y^2$")

ax[0, 1].plot(N, accelerationExpression(pathx, pathy, pathz))


reversePath = integrate.odeint(vectorPrime, [x0, y0, z0], reverseN)
reversePathTranspose = reversePath.transpose()

reversePathx = reversePathTranspose[0]
reversePathy = reversePathTranspose[1]
reversePathz = reversePathTranspose[2]

ax[0, 0].plot(reversePathx, reversePathy, reversePathz, linewidth=0.75)

ax[1, 0].plot(reverseN, gamma_phi(reversePathx, reversePathy))

ax[1, 1].plot(reverseN, reversePathz**2,                           "r")
ax[1, 1].plot(reverseN, 1 - reversePathx**2 - reversePathy**2 - reversePathz**2, "g")
ax[1, 1].plot(reverseN, reversePathx**2 + reversePathy**2,                "b")

ax[0, 1].plot(reverseN, accelerationExpression(reversePathx, reversePathy, reversePathz))

fixedPoints = fixedPoints.transpose()

ax[0, 0].scatter(fixedPoints[0], fixedPoints[1], fixedPoints[2], depthshade=False)

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
ax[1, 0].set_ylabel("$\gamma_\phi$")

ax[1, 1].set_xlabel("N")
ax[1, 1].set_ylabel("Density Parameters")
ax[1, 1].legend()
plt.show()
