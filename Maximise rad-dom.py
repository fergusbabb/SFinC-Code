import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt

def vectorPrime(state, N, lam):
    
    x = state[0]
    y = state[1]
    z = state[2]

    HdotSection = 1 + x**2 - y**2 + (1/3)*z**2

    xprime = (-3*x) + (lam * np.sqrt(1.5) * y**2)  + (1.5*x * HdotSection)
    yprime =         (-lam * np.sqrt(1.5) * x * y) + (1.5*y * HdotSection)
    zprime = (-2*z)                                + (1.5*z * HdotSection)

    return [xprime, yprime, zprime]

lam = np.linspace(0.3, 0.5, 400)
peakOmega_r = []

Omega_phi = 0.73
gamma_phi = 0.013

Omega_r = 4.984e-5

xSquared = gamma_phi / (2 * Omega_phi)
ySquared = Omega_phi - xSquared

state_0 = [np.sqrt(xSquared),
           np.sqrt(ySquared),
           np.sqrt(Omega_r)]

N = np.linspace(0, -15, 101) 

for lamValue in lam:
    path = integrate.odeint(vectorPrime, state_0, N, args=(lamValue,))
    pathz = path.transpose()[2]
    peakz = np.max(pathz)
    peakOmega_r.append(peakz**2)

    plt.scatter(lamValue, peakz**2)
##    plt.plot(N, pathz**2)

print(np.max(peakOmega_r))
plt.plot(lam, peakOmega_r)
plt.ylim([0, 1])
plt.show()
