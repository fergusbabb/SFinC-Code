import numpy as np

def gamma_phi(x, y):
    return (2*x**2) / (x**2 + y**2)

x = 0.002
y1 = 0.001
y2 = 0.001

y = np.sqrt(y1**2 + y2**2)

print(gamma_phi(x,y))