from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
from matplotlib.colors import Normalize
import matplotlib.cm as cm

import ctypes
ctypes.windll.shcore.SetProcessDpiAwareness(1)

#Personal plotting preferences
plt.rcParams.update({"text.usetex": True, "font.family": "serif",
                     "font.serif": ["Computer Modern Serif"]})
plt.rc('axes', labelsize=16, titlesize=18)
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

#_________________________Set up main window___________________________________

window = tk.Tk()
window.title('Autonomous Systems')
window.geometry('900x900')



fig = Figure(figsize=(6, 3))
axis_dims = [.1,.2,.6,.7]
ax = fig.add_axes(axis_dims)
ax.set_aspect('equal')

#Colour bar axes
cbar_ax_dims = [.75,.575,.015,.35]
cbar_ax = fig.add_axes(cbar_ax_dims)

#Bounding Circle
theta = np.linspace(0, np.pi, 150)
a = np.cos(theta)
b = np.sin(theta)
ax.plot(a, b, 'r', linewidth=1, label='Bounding Circle')
ax.plot([-1,1], [0,0], 'r', linewidth=1)


gam_0 = 1
lam_0 = 3
pathnum = 50
N = np.linspace(0, 8, 1000) 
xinit = np.linspace(-0.99, 0.99, pathnum)


ax.set_xlabel('$x$', x=1.02)
ax.set_ylabel('$y$', rotation = 0, y=1.02)

plot_label = tk.Label(window, text = "Master Plot", width = 10, 
                       height = 2).pack(anchor='n')

canvas = FigureCanvasTkAgg(fig, window) 
canvas.draw()
canvas.get_tk_widget().pack(anchor = 'center', pady=0.2)



# Define the system of ODEs
def ODEs(coords, t, lam, gam):
    # x[0] corresponds to the variable x
    # x[1] corresponds to the variable y
    x, y = coords[0], coords[1]

    # Return the derivatives [x', y']
    dx_dN = -3*x + lam*(np.sqrt(3/2))*y**2 + (3/2)*x*(2*x**2+gam*(1-x**2-y**2))
    dy_dN = -lam*(np.sqrt(3/2))*y*x + (3/2)*y*(2*x**2+gam*(1-x**2-y**2))
    return [dx_dN, dy_dN]


def update_plot(event):
    lam = lambda_slide.get()
    gam = gamma_slide.get()

    for i in range(pathnum):
        initial_conditions = [xinit[i], 0.01]  # Initial values for x and y
    
        # Solve the system of ODEs using odeint
        solution = odeint(ODEs, initial_conditions, N, (lam, gam))
        main_tracks[i].set_data(solution[:,0], solution[:,1])

        #Update the quiver vectors
        quiver_vectors = np.array([ODEs([pt[0], pt[1]], N, lam, gam)
                                    for pt in points])
        u, v = quiver_vectors[:, 0], quiver_vectors[:, 1]
        magnitude = np.sqrt(u**2 + v**2).flatten()
        u = u/magnitude
        v = v/magnitude

        cmap = cm.spring

        global quiver
        quiver.remove()
        quiver = ax.quiver(x_ins, y_ins, u, v,
                    color=cmap(norm(magnitude)), norm=norm,
                    cmap=cmap)
    fig.canvas.draw()


def fibonacci_semicircle(max_radius, num_points):
    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # Golden angle in radians

    for i in range(num_points):
        radius = max_radius * np.sqrt(i / num_points)  # Scale radius non-linearly
        theta = phi * i  # Increment by the golden angle

        x = np.cos(theta) * radius
        y = np.sin(theta) * radius

        # Only keep points in the upper half (y >= 0)
        if y >= 0.025:
            points.append([x, y])

    return np.array(points)


# Parameters
radius = 0.95
num_points = 400  # Number of points on the semicircle

# Generate points on a semicircle using Fibonacci-like spacing
points = fibonacci_semicircle(radius, num_points)
points = np.vstack((points,[0,0]))

# Compute vectors at each point using the example vector field function
quiver_vectors = np.array([ODEs([pt[0], pt[1]], N, lam_0, gam_0)
                            for pt in points])
u, v = quiver_vectors[:, 0], quiver_vectors[:, 1]
x_ins, y_ins = points[:, 0], points[:, 1]

#Scale for the colourbar
magnitude = np.sqrt(u**2 + v**2).flatten()
norm = Normalize()
u = u/magnitude
v = v/magnitude

cmap = cm.spring

#Plot the Quivers
quiver = ax.quiver(x_ins, y_ins, u, v,
                    color=cmap(norm(magnitude)), norm=norm,
                    cmap=cmap, scale=25, headwidth = 3)

#Plot the colourbar
cbar = plt.colorbar(quiver, cax=cbar_ax, orientation='vertical')
cbar.set_label('Magnitude of derivatives')



lambda_slide_label = tk.Label(window, text = 'r$\lambda$ value', width = 15, 
                       height = 2).pack(anchor = 'center', pady=0.1)
lambda_slide = tk.Scale(window, from_ = 0, to = 5, orient = 'horizontal',
                       width = 20, length = 250, resolution=0.01)
lambda_slide.bind("<ButtonRelease-1>", update_plot)
lambda_slide.set(lam_0)
lambda_slide.pack(anchor = 'center', pady=0.1)

gamma_slide_label = tk.Label(window, text = 'r$\gamma$ value', width = 15, 
                       height = 2).pack(anchor = 'center', pady=0.1)
gamma_slide = tk.Scale(window, from_ = 0, to = 5, orient = 'horizontal',
                       width = 20, length = 250, resolution=0.01)
gamma_slide.bind("<ButtonRelease-1>", update_plot)
gamma_slide.set(gam_0)
gamma_slide.pack(anchor = 'center', pady=0.1)


main_tracks = []
for i in range(pathnum):
    initial_conditions = [xinit[i], 0.01]  # Initial values for x and y
    
    # Solve the system of ODEs using odeint
    solution = odeint(ODEs, initial_conditions, N, (lam_0, gam_0))
    track_i = ax.plot(solution[:,0], solution[:,1], 'k', linewidth=.5)[0]
    
    main_tracks.append(track_i)
    track_i.set_visible(True)






window.mainloop()
