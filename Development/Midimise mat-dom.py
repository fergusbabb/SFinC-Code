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

lam = np.linspace(0, 10, 101)
lambdaConsistent = np.array([])

Omega_phi = 0.73
gamma_phi = 0.013

Omega_r = 4.984e-5
peakOmega_mConsistent = np.array([])

xSquared = gamma_phi / (2 * Omega_phi)
ySquared = Omega_phi - xSquared

state_0 = [np.sqrt(xSquared),
           np.sqrt(ySquared),
           np.sqrt(Omega_r)]

N = np.linspace(0, -25, 301) 

def simulate(lam):
    peakOmega_m = []
    for lamValue in lam:
        path = integrate.odeint(vectorPrime, state_0, N, args=(lamValue,))
        pathx = path.transpose()[0]
        pathy = path.transpose()[1]
        pathz = path.transpose()[2]

        pathMatter = 1 - pathx**2 - pathy**2 - pathz**2        
        peakMatter = np.max(pathMatter)
        peakOmega_m.append(peakMatter)

##        plt.scatter(lamValue, peakz**2)
        
    return peakOmega_m

while lam[-1] - lam[0] > 1e-10:
    peaks = simulate(lam)
    peakOmega_mConsistent = np.append(peakOmega_mConsistent, peaks)
    maxIndex = np.argmax(peaks)
    print("Index: " + str(maxIndex) + "\t Omega_r: " + str(peaks[maxIndex]) + "\t lambda: " + str(lam[maxIndex]))

    lambdaConsistent = np.append(lambdaConsistent, lam)
    lam = np.linspace(lam[maxIndex - 1], lam[maxIndex + 1], 11)

order = np.argsort(lambdaConsistent)

plt.plot(lambdaConsistent[order], peakOmega_mConsistent[order])
plt.ylim([0, 1])
plt.show()
