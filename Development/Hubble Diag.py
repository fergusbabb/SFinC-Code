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

def d_L_Integrand(currentTotal, z, Omega_m0, Omega_DE):#, pathx, pathy, pathz):
##    Omega_DE = 0.7
##    Omega_m0 = 0.3
    return 1/np.sqrt(1 - Omega_DE + Omega_m0*(1+z)**3)

lam = 0.38583349
lambdaConsistent = np.array([])

Omega_phi = 0.73
gamma_phi = 0.013

Omega_r = 4.984e-5
peakOmega_rConsistent = np.array([])

xSquared = gamma_phi / (2 * Omega_phi)
ySquared = Omega_phi - xSquared

Omega_m0 = 1 - xSquared - ySquared - Omega_r

state_0 = [np.sqrt(xSquared),
           np.sqrt(ySquared),
           np.sqrt(Omega_r)]

N = np.linspace(0, -1.4, 31)

z = np.exp(-N) - 1
c = 3e8
V = z * c

##path = integrate.odeint(vectorPrime, state_0, N, args=(lam,))
##pathTranspose = path.transpose()
##
##pathx = pathTranspose[0]
##pathy = pathTranspose[1]
##pathz = pathTranspose[2]
##
##plt.plot(N, pathz**2,                           "r", label = "$\Omega_r = z^2$")
##plt.plot(N, 1 - pathx**2 - pathy**2 - pathz**2, "g", label = "$\Omega_m = 1 - x^2 - y^2 - z^2$")
##plt.plot(N, pathx**2 + pathy**2,                "b", label = "$\Omega_\phi = x^2 + y^2$")
##plt.show()


##plt.plot(z, V)
##plt.show()

for Omega_DE in [0, 0.3, 0.7, 1]:
    plt.plot(z, (1 + z) * integrate.odeint(d_L_Integrand, 0, z, args=(Omega_m0, Omega_DE)).transpose()[0])
plt.show()
