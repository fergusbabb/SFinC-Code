import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import (odeint, quad)

#_____________________________Define ODE for 2 potentials_______________________

def ODEs(state, N, lam1, lam2):
    # state[0]   corresponds to the variable x
    # state[1]   corresponds to the variable y1
    # state[2]   corresponds to the variable y2
    # state[3]   corresponds to the variable z

    x = state[0]
    y1= state[1]
    y2= state[2]
    z = state[3]

    HdotSection = 1 + x**2 - y1**2 - y2**2 + (1/3)*z**2

    xprime = (-3*x) + (np.sqrt(1.5) * (lam1*y1**2 + lam2*y2**2)) + (1.5*x * HdotSection)
    y1prime=        - (np.sqrt(1.5) *  lam1 * x *y1)             + (1.5*y1* HdotSection)
    y2prime=        - (np.sqrt(1.5) *  lam2 * x *y2)             + (1.5*y2* HdotSection)
    zprime = (-2*z)                                              + (1.5*z * HdotSection)

    return [xprime, y1prime, y2prime, zprime]

N = np.linspace(0, -25, 301)

Omega_phi = 0.73
gamma_phi = 0.013

Omega_r = 5e-5
peakOmega_rConsistent = np.array([])

xSquared = gamma_phi / (2 * Omega_phi)
ySquared = Omega_phi - xSquared

state_0 = np.array([np.sqrt(xSquared),
           np.sqrt(ySquared/2),
           np.sqrt(ySquared/2),
           np.sqrt(Omega_r)])


def simulate(lam1, lam2, N, state_0):
    grid = np.zeros((len(lam1), len(lam2)))
    for i, lam1Value in enumerate(lam1[0]):        
        print(i)
        for j, lam2Value in enumerate(lam2[:,0]):
            path = odeint(ODEs, state_0, N, args=(lam1Value, lam2Value))
            pathz = path.transpose()[3]

            grid[i, j] = np.max(pathz)
##            print(str(i) + ": " + str(lam1Value) + ", " + str(j) + ": " + str(lam2Value) + "\t: " + str(grid[i, j]))

    return grid

lambda1 = np.linspace(-2, 2, 201)
lambda2 = np.linspace(-2, 2, 201)

lam1, lam2 = np.meshgrid(lambda1, lambda2)


radiation = simulate(lam1, lam2, N, state_0)

contour = plt.contour(lam1, lam2, radiation)
plt.clabel(contour, inline=True, fontsize=10)
plt.show()
