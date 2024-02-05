from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)

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
window.geometry('900x1200')



fig = Figure(figsize=(6, 3))
axis_dims = [.15,.20,.7,.7]
ax = fig.add_axes(axis_dims)
ax.set_aspect('equal')

fig2 = Figure(figsize=(6, 3))
axis_dims2 = [.15,.20,.7,.7]
ax2 = fig2.add_axes(axis_dims2)
#ax2.set_aspect('equal')

#Bounding Circle
theta = np.linspace(0, np.pi, 150)
a = np.cos(theta)
b = np.sin(theta)
ax.plot(a, b, 'r', linewidth=1, label='Bounding Circle')
ax.plot([-1,1], [0,0], 'r', linewidth=1)


gam_0 = 1
lam_0 = 3
pathnum = 50



ax.set_xlabel('$x$', x=1.02)
ax.set_ylabel('$y$', rotation = 0, y=1.02)

plot_label = tk.Label(window, text = "Master Plot", width = 10, 
                       height = 2).pack(anchor='n')

canvas = FigureCanvasTkAgg(fig, window) 
canvas.draw()
canvas.get_tk_widget().pack(anchor = 'center', pady=0.2)

canvas2 = FigureCanvasTkAgg(fig2, window)
canvas2.draw()
canvas2.get_tk_widget().pack(anchor = 'center', pady=0.2)




# Define the system of ODEs
def ODEs(coords, t):
    # x[0] corresponds to the variable x
    # x[1] corresponds to the variable y
    x, y = coords[0], coords[1]

    lam = lambda_slide.get()
    gam = gamma_slide.get()

    # Return the derivatives [x', y']
    dx_dN = -3*x + lam*(np.sqrt(3/2))*y**2 + (3/2)*x*(2*x**2+gam*(1-x**2-y**2))
    dy_dN = -lam*(np.sqrt(3/2))*y*x + (3/2)*y*(2*x**2+gam*(1-x**2-y**2))
    return [dx_dN, dy_dN]


def update_plot(event):
    for i in range(pathnum):
        initial_conditions = [xinit[i], 0.01]  # Initial values for x and y
    
        # Solve the system of ODEs using odeint
        solution = odeint(ODEs, initial_conditions, n)
        main_tracks[i].set_data(solution[:,0], solution[:,1])
    lam = lambda_slide.get()
    gam = gamma_slide.get()
    clickpoint.set_data(lam**2,gam)
    fig.canvas.draw()


ax2.set(xlim=[0,12],ylim=[0,2])
ax2.plot([0,6],[12,2],'k', linestyle = '-')
ax2.plot([0,12],[0,0],'k')
ax2.plot([0,0],[0,2],'k')
ax2.plot([6,6],[0,2],'k')
ax2.plot([0,6],[0,2],'k')
ax2.plot([2,2],[0,2],'k', linestyle = ':')
clickpoint, = ax2.plot(lam_0**2,gam_0,'r*')

def regions_plot(event):
    ix, iy = event.xdata, event.ydata
    lambda_slide.set(np.sqrt(ix))
    gamma_slide.set(iy)
    clickpoint.set_data(ix,iy)
    update_plot(event)
    fig2.canvas.draw()

cid = fig2.canvas.mpl_connect('button_press_event', regions_plot)
lambda_slide_label = tk.Label(window, text = 'r$\lambda$ value', width = 15, 
                       height = 2).pack(anchor = 'center', pady=0.1)
lambda_slide = tk.Scale(window, from_ = 0, to = np.sqrt(12), orient = 'horizontal',
                       width = 20, length = 250, resolution=0.01)
lambda_slide.bind("<ButtonRelease-1>", update_plot)
lambda_slide.set(lam_0)
lambda_slide.pack(anchor = 'center', pady=0.1)

gamma_slide_label = tk.Label(window, text = 'r$\gamma$ value', width = 15, 
                       height = 2).pack(anchor = 'center', pady=0.1)
gamma_slide = tk.Scale(window, from_ = 0, to = 2, orient = 'horizontal',
                       width = 20, length = 250, resolution=0.01)
gamma_slide.bind("<ButtonRelease-1>", update_plot)
gamma_slide.set(gam_0)
gamma_slide.pack(anchor = 'center', pady=0.1)


n = np.linspace(0, 15, 1000) 
xinit = np.linspace(-0.99, 0.99, pathnum)

main_tracks = []
for i in range(pathnum):
    initial_conditions = [xinit[i], 0.01]  # Initial values for x and y
    
    # Solve the system of ODEs using odeint
    solution = odeint(ODEs, initial_conditions, n)
    track_i = ax.plot(solution[:,0], solution[:,1], 'k', linewidth=.5)[0]
    
    main_tracks.append(track_i)
    track_i.set_visible(True)






window.mainloop()
