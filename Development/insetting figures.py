import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Example data
N = np.linspace(0, 10, 100)
Densities = np.exp(N)

# Create main figure and axis
fig, rho_ax = plt.subplots()

# Plot example data on main axis
rho_ax.plot(N, Densities)
rho_ax.set(ylim=[1e0, 1e30],
           yticks=[1e0, 1e5, 1e10, 1e15, 1e20, 1e25, 1e30],
           yticklabels=['$10^{0}$', '$10^{5}$', '$10^{10}$', '$10^{15}$', '$10^{20}$', '$10^{25}$', '$10^{30}$'],
           xticks=[0, 2, 4, 6, 8, 10])

# Add inset axes (zoomed in section)
axins = inset_axes(rho_ax, width="40%", height="40%", loc='upper right', borderpad=2)

# Define the region of interest
x1, x2, y1, y2 = 2, 4, 1e1, 1e20  # Limits for the zoomed in section
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

# Plot the same data on the inset axes
axins.plot(N, Densities)

# Optional: Mark the region of the main plot being zoomed in
rho_ax.indicate_inset_zoom(axins, edgecolor="black")

plt.show()