import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
n_points = 10

# Define your vector field function (example)
def vector_field_function(x, y, z):
    # Example function: vector field points in the direction of (x, y, z)
    return x, y, z


'''
theta = np.linspace(0, np.pi, n_points)
phi = np.linspace(0, np.pi / 2, n_points)
r = np.linspace(0, 1, n_points)  # Radius for points inside the sphere

# Create a meshgrid of spherical coordinates
theta_grid, phi_grid, r_grid = np.meshgrid(theta, phi, r, indexing='ij')

# Convert spherical coordinates to Cartesian coordinates
z = r_grid * np.sin(theta_grid) * np.cos(phi_grid)
y = r_grid * np.sin(theta_grid) * np.sin(phi_grid)
x = r_grid * np.cos(theta_grid)

# Filter out points outside the sphere
inside_sphere = x**2 + y**2 + z**2 <= 1
x_inside = x[inside_sphere]
y_inside = y[inside_sphere]
z_inside = z[inside_sphere]

# Evaluate the vector field at each point inside the sphere
u, v, w = vector_field_function(x_inside, y_inside, z_inside)
ax.quiver(x_inside, y_inside, z_inside, u, v, w, length=0.05, normalize=True)

'''


# Generate a regular grid of points in a square
x = np.linspace(-1, 1, n_points)
y = np.linspace(0, 1, n_points)
z = np.linspace(0, 1, n_points)

# Create a meshgrid of points
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Flatten the arrays
X = X.flatten()
Y = Y.flatten()
Z = Z.flatten()

# Filter out points outside the sphere
inside_sphere = X**2 + Y**2 + Z**2 <= 1
X_inside = X[inside_sphere]
Y_inside = Y[inside_sphere]
Z_inside = Z[inside_sphere]

# Evaluate the vector field at each point inside the sphere
u, v, w = vector_field_function(X_inside, Y_inside, Z_inside)
ax.quiver(X_inside, Y_inside, Z_inside, u, v, w, length=0.05, normalize=True)




ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_aspect('equal')

theta = np.linspace(0, np.pi, 150)
ax.plot(np.cos(theta), np.sin(theta),
         'r', linewidth=1, label='Bounding Circle')
ax.plot(np.cos(theta), np.zeros(150), np.sin(theta),
         'r', linewidth=1)
ax.plot(np.zeros(75),np.cos(theta[0:75]),np.sin(theta[0:75]),
         'r', linewidth=1)
ax.plot([-1,0,1], [0,0,0], [0,0,0], 'r', linewidth=1)
ax.plot([0,0,0], [1,0,0], [0,0,1], 'r', linewidth=1)


plt.show()