import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize
import matplotlib.cm as cm

# Define your ODE system for two fluids
def ODEs(state, t, lam):
    x, y, z = state
    xprime = (-3*x) + (lam * np.sqrt(1.5) * y**2) + (1.5 * x * (1 + x**2 - y**2 + (1/3)*z**2))
    yprime = (-lam * np.sqrt(1.5) * x * y) + (1.5 * y * (1 + x**2 - y**2 + (1/3)*z**2))
    zprime = (-2*z) + (1.5 * z * (1 + x**2 - y**2 + (1/3)*z**2))
    return [xprime, yprime, zprime]

n_points = 10
x = np.linspace(-1, 1, n_points)
y = np.linspace(0, 1, n_points)
z = np.linspace(0, 1, n_points)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
X = X.flatten()
Y = Y.flatten()
Z = Z.flatten()

# Filter points inside the sphere
inside_sphere = X**2 + Y**2 + Z**2 <= 1
X_inside = X[inside_sphere]
Y_inside = Y[inside_sphere]
Z_inside = Z[inside_sphere]

# Vector field at the points
state_quiver = np.vstack((X_inside, Y_inside, Z_inside)).T
t = 0  # If your system does not explicitly depend on time, this can be a dummy value.
lam = 1  # Example lambda value

# Use the ODE function
results = np.array([ODEs(state, t, lam) for state in state_quiver])
u, v, w = results[:, 0], results[:, 1], results[:, 2]

# Color by magnitude
magnitude = np.sqrt(u**2 + v**2 + w**2)
norm = Normalize(vmin=magnitude.min(), vmax=magnitude.max())
cmap = cm.plasma

# Plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
quiver = ax.quiver(X_inside, Y_inside, Z_inside, u, v, w, length=0.1, cmap=cmap, color=cmap(norm(magnitude)), norm=norm)
cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical')
cbar.set_label('Magnitude of derivatives')
plt.show()