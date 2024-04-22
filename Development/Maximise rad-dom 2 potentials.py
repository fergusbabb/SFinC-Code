import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt

def vectorPrime(state, N, lam):
    x = state[0]
    y = state[1:-1]
    z = state[-1]

    HdotSection = 1 + x**2 - np.sum(y**2) + (1/3)*z**2

    xprime = (-3*x) + (np.sqrt(1.5) *  np.sum(lam * y**2)) + (1.5*x * HdotSection)
    yprime =          (np.sqrt(1.5) *        -lam * x * y) + (1.5*y * HdotSection)
    zprime = (-2*z)                                        + (1.5*z * HdotSection)

    return np.append(xprime, np.append(yprime, zprime))


Omega_phi = 0.73
gamma_phi = 0.013

Omega_r = 4.984e-5
peakOmega_rConsistent = np.array([])

xSquared = gamma_phi / (2 * Omega_phi)
ySquared = Omega_phi - xSquared

state_0 = [np.sqrt(xSquared),
           np.sqrt(ySquared * (99999/100000)),
           np.sqrt(ySquared * (1/100000)),
           np.sqrt(Omega_r)]

N = np.linspace(0, -25, 301) 

def simulate(lam):
##    print(lam)
    peakOmega_r = []
    for lamValue in lam:
        path = integrate.odeint(vectorPrime, state_0, N, args=(lamValue,))
        pathx = path.transpose()[0]
        pathy = path.transpose()[1]
        pathz = path.transpose()[2]

        pathMatter = 1 - pathx**2 - pathy**2 - pathz**2        
        peakMatter = np.max(pathMatter)
        peakOmega_r.append(peakMatter)

##        plt.scatter(lamValue, peakz**2)
        
    return peakOmega_r

for lam2 in [0]:#np.linspace(0, 0.1, 11):
    lam = np.linspace((0, lam2), (1, lam2), 11)
    lambdaConsistent = np.array([])

    while lam[-1][0] - lam[0][0] > 1e-10:
        peaks = simulate(lam)
        peakOmega_rConsistent = np.append(peakOmega_rConsistent, peaks)
        maxIndex = np.argmax(peaks)
        print("Index: " + str(maxIndex) + "\t Omega_r: " + str(peaks[maxIndex]) + "\t lambda: " + str(lam[maxIndex]))

        lambdaConsistent = np.append(lambdaConsistent, lam[0])
        lam = np.linspace((lam[0][maxIndex - 1], lam2), (lam[0][maxIndex + 1], lam2), 11)

    order = np.argsort(lambdaConsistent)

    plt.plot(lambdaConsistent[order], peakOmega_rConsistent[order])
    print()
    
plt.ylim([0, 1])
plt.show()
