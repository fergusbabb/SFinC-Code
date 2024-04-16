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

def d_L_Integrand(currentTotal, z, zaxis, Omega_m0, Omega_r0, Omega_phi_0, path_gamma_phi):
    index = np.abs(zaxis-z).argmin()
##    print(index)
    return 1/np.sqrt(Omega_m0*(1+z)**3 +
                     Omega_r0*(1+z)**4 + 
                     Omega_phi_0*(1+z)**(3 * path_gamma_phi[index]))

lam = 0.38583349
lambdaConsistent = np.array([])

Omega_phi_0 = 0.73
gamma_phi = 0.013

Omega_r0 = 4.984e-5
peakOmega_rConsistent = np.array([])

xSquared = gamma_phi / (2 * Omega_phi_0)
ySquared = Omega_phi_0 - xSquared

Omega_m0 = 1 - xSquared - ySquared - Omega_r0

state_0 = [np.sqrt(xSquared),
           np.sqrt(ySquared),
           np.sqrt(Omega_r0)]

N = np.linspace(0, -1, 31)

z = np.exp(-N) - 1
c = 3e5     # Given in km/s
V = z * c   # ""
h = 0.738
H_0 = 100*h # km/s/Mpc

path = integrate.odeint(vectorPrime, state_0, N, args=(lam,))
pathTranspose = path.transpose()

pathx = pathTranspose[0]
pathy = pathTranspose[1]
pathz = pathTranspose[2]

##plt.plot(N, pathz**2,                           "r", label = "$\Omega_r = z^2$")
##plt.plot(N, 1 - pathx**2 - pathy**2 - pathz**2, "g", label = "$\Omega_m = 1 - x^2 - y^2 - z^2$")
##plt.plot(N, pathx**2 + pathy**2,                "b", label = "$\Omega_\phi = x^2 + y^2$")
##plt.show()


##plt.plot(z, V)
##plt.show()

path_gamma_phi = (2 * pathx**2) / (pathx**2 + pathy**2)

##for Omega_phi_0 in [0, 0.3, 0.7, 1]:
d_L = (c/H_0) * (1 + z) * integrate.odeint(d_L_Integrand, 0, z, args=(z, Omega_m0, Omega_r0, Omega_phi_0, path_gamma_phi)).transpose()[0]

plt.plot(d_L, V)
plt.show()
