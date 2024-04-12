import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt

lam = np.array([1, 3])

## Initial conditions
x0 = -0.01
y0 = 0.01 / np.sqrt(len(lam))
z0 = 0.9

start = np.append(x0, np.append([y0 for lambda_n in lam], z0))

print(np.sum(start**2))


fig, ax = plt.subplots(2, 2)

ax[0, 0].set_axis_off()
ax[0, 0] = fig.add_subplot(2, 2, 1, projection="3d")


def vectorPrime(state, N):
    global lam
    
    x = state[0]
    y = state[1:-1]
    z = state[-1]

    HdotSection = 1 + x**2 - np.sum(y**2) + (1/3)*z**2

    xprime = (-3*x) + (np.sqrt(1.5) *  np.sum(lam * y**2)) + (1.5*x * HdotSection)
    yprime =          (np.sqrt(1.5) *        -lam * x * y) + (1.5*y * HdotSection)
    zprime = (-2*z)                                        + (1.5*z * HdotSection)

    return np.append(xprime, np.append(yprime, zprime))

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

N = np.linspace(0, 20, 401)
reverseN = np.linspace(0, -8, 101)
a = np.exp(N)

ax[0, 1].plot([reverseN[-1], N[-1]], [0, 0], "k--")

ax[1, 0].plot([reverseN[-1], N[-1]], [4/3, 4/3], "k--", linewidth = 0.5)
ax[1, 0].plot([reverseN[-1], N[-1]], [1, 1], "k--", linewidth = 0.5)


## Initial conditions

##if lam[0] > 2:
##    y0 = 2 / (lam * np.sqrt(3))#0.2
##    x0 = 2 * np.sqrt(2/3) / lam #np.sqrt(2) * y0
##    z0 = np.sqrt(1 - (4/lam**2)) - 0.005
##else:

fixedPoints = np.array([
    [0, 0, 0],
    [0, 0, 1],
    [1, 0, 0],
    [-1, 0, 0]
    ])

##if lam**2 < 6:
##    fixedPoints = np.append(fixedPoints, [[lam/np.sqrt(6), np.sqrt(1 - (lam**2)/6), 0]], axis=0)
##
##if lam**2 > 3:
##    fixedPoints = np.append(fixedPoints, [[np.sqrt(3/2)/lam, np.sqrt(3/2)/lam, 0]], axis=0)
##
##if lam**2 > 8/3:
##    fixedPoints = np.append(fixedPoints, [[2 * np.sqrt(2/3) / lam, 2 / (lam * np.sqrt(3)), np.sqrt(1 - (4/lam**2))]], axis=0)

## FORWARD PATH

path = integrate.odeint(vectorPrime, start, N)
pathTranspose = path.transpose()

pathx   = pathTranspose[0]
pathy_n = pathTranspose[1:-1]
pathz   = pathTranspose[-1]

pathy_total = np.sqrt(np.sum(pathy_n**2, axis=0))

## PLOT FORWARD

ax[0, 0].plot(pathx, pathy_total, pathz, linewidth=0.75)
ax[1, 0].plot(N, gamma_phi(pathx, pathy_total))

for pathy in pathy_n:
    ax[1, 1].plot(N, pathy**2, "b--", linewidth=0.5)

ax[1, 1].plot(N, pathz**2,                                 "r", label = "$\Omega_r = z^2$")
ax[1, 1].plot(N, 1 - pathx**2 - pathy_total**2 - pathz**2, "g", label = "$\Omega_m = 1 - x^2 - y^2 - z^2$")
ax[1, 1].plot(N, pathx**2 + pathy_total**2,                "b", label = "$\Omega_\phi = x^2 + y^2$")

ax[0, 1].plot(N, accelerationExpression(pathx, pathy_total, pathz))

## REVERSE PATH

reversePath = integrate.odeint(vectorPrime, start, reverseN)
reversePathTranspose = reversePath.transpose()

reversePathx = reversePathTranspose[0]
reversePathy_n = reversePathTranspose[1:-1]
reversePathz = reversePathTranspose[-1]

reversePathy_total = np.sqrt(np.sum(reversePathy_n**2, axis=0))

## PLOT REVERSE

ax[0, 0].plot(reversePathx, reversePathy_total, reversePathz, linewidth=0.75)
ax[1, 0].plot(reverseN, gamma_phi(reversePathx, reversePathy_total))

for pathy in reversePathy_n:
    ax[1, 1].plot(reverseN, pathy**2, "b--", linewidth=0.5)

ax[1, 1].plot(reverseN, reversePathz**2,                                               "r")
ax[1, 1].plot(reverseN, 1 - reversePathx**2 - reversePathy_total**2 - reversePathz**2, "g")
ax[1, 1].plot(reverseN, reversePathx**2 + reversePathy_total**2,                       "b")

ax[0, 1].plot(reverseN, accelerationExpression(reversePathx, reversePathy_total, reversePathz))

##

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
