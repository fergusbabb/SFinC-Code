from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

gam = 1
lam = 1
pathnum = 25


# Define the system of ODEs
def ODEs(coords, t):
    # x[0] corresponds to the variable x
    # x[1] corresponds to the variable y
    x, y = coords[0], coords[1]
    
    # Return the derivatives [x', y']
    dx_dN = -3*x + lam*(np.sqrt(3/2))*y**2 + (3/2)*x*(2*x**2+gam*(1-x**2-y**2))
    dy_dN = -lam*(np.sqrt(3/2))*y*x + (3/2)*y*(2*x**2+gam*(1-x**2-y**2))
    return [dx_dN, dy_dN]

n = np.linspace(0, 7, 1000)  # You can adjust the time points as needed
# Set initial conditions

xinit = np.linspace(-0.99,0.99,pathnum)
for i in range(pathnum):
    initial_conditions = [xinit[i], 0.01]  # Initial values for x and y
    
    # Solve the system of ODEs using odeint
    solution = odeint(ODEs, initial_conditions, n)
    plt.plot(solution[:,0], solution[:,1], 'k', linewidth=.5)

theta = np.linspace(0, np.pi, 150)
a = np.cos(theta)
b = np.sin(theta)
plt.plot(a, b, 'r', linewidth=1, label='Bounding Circle')

plt.show()