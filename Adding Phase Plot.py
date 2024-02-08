#_______________________________Imports________________________________________
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)

#Tkinter for widgets and interactivity. Scipy for ode solving
import tkinter as tk
from scipy.integrate import odeint

#This import stops automatically sets up the window so its reasonable
import ctypes 
ctypes.windll.shcore.SetProcessDpiAwareness(1) 


#Personal plotting preferences
'''plt.rcParams.update({"text.usetex": True, "font.family": "serif",
                     "font.serif": ["Computer Modern Serif"]})'''
plt.rc('axes', labelsize=12, titlesize=15)
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

#_________________________Set up main window___________________________________
#Standard tkinter window set up
window = tk.Tk()
window.title('Autonomous Systems')
window.geometry('900x1200')

#Tracks plot figure
fig = Figure(figsize=(6, 3))
axis_dims = [.15,.20,.7,.7]
ax = fig.add_axes(axis_dims)
ax.set_aspect('equal')

#Interactivity phase plot figure
fig2 = Figure(figsize=(6, 3))
axis_dims2 = [.15,.20,.7,.7]
ax2 = fig2.add_axes(axis_dims2)
#ax2.set_aspect('equal')

#Bounding Circle
theta = np.linspace(0, np.pi, 150)
ax.plot(np.cos(theta), np.sin(theta),
         'r', linewidth=1, label='Bounding Circle')
ax.plot([-1,1], [0,0], 'r', linewidth=1)


#Initial values
gam_0 = 1 
lam_0 = 3
pathnum = 10


#_________________________Initialise plots___________________________________
#Tkinter works as:
#1) Main Window 
#2) Set up canvases inside window 
#3) Set up Figures attached to canvases 
#4) Add axes to figures as usual

ax.set_xlabel('$x$', x=1.02)
ax.set_ylabel('$y$', rotation = 0, y=1.02)

plot_label = tk.Label(window, text = "Master Plot", width = 10, 
                       height = 2).pack(anchor='n')

canvas = FigureCanvasTkAgg(fig, window) 
#Canvas is where figure is placed to window
canvas.draw() #Show canvas (ie show figure)
canvas.get_tk_widget().pack(anchor = 'center', pady=0.2)

#Interactive plot, plot lines as in P40
ax2.set(xlim=[0,12],ylim=[0,2])
ax2.set_xlabel('$\lambda^2$', x=1.02)
ax2.set_ylabel('$\gamma$', rotation = 0, y=1.02)
ax2.plot([0,6],[12,2],'k', linestyle = '-')
ax2.plot([0,12],[0,0],'k')
ax2.plot([0,0],[0,2],'k')
ax2.plot([6,6],[0,2],'k')
ax2.plot([0,6],[0,2],'k')
ax2.plot([2,2],[0,2],'k', linestyle = ':')


#As before for interactive phase plot
canvas2 = FigureCanvasTkAgg(fig2, window)
canvas2.draw()
canvas2.get_tk_widget().pack(anchor = 'center', pady=0.2)




#_____________________________Define ODE for 1 fluid_____________________
def ODEs(coords, t, lam, gam):
    # coords[0] corresponds to the variable x
    # coords[1] corresponds to the variable y
    x, y = coords[0], coords[1]

    A = 2*x**2+gam*(1-x**2-y**2)
    a = np.sqrt(3/2)

    #Return the derivatives [x', y']
    dx_dN = -3*x + lam*a*y**2 + a**2*x*A
    dy_dN =      - lam*a*y*x  + a**2*y*A

    return [dx_dN, dy_dN]

#________________________Define Acceleration Expression______________________
def acceleration(x, y):
    return -3*x**2 - 1.5*gam*(1 - x**2 - y**2) + 1

#_______________________________Update track plots___________________________
def update_plot(event):
    #Before update_plot is called lam and gam sliders are updated
    #Sliders will always have current lambda/ gamma values
    lam = lambda_slide.get()
    gam = gamma_slide.get()

    #Update cursor star point
    clickpoint.set_data(lam**2, gam)

    #Plot all paths with updated values
    for i in range(pathnum):
        initial_conditions = [xinit[i], 0.01]  #Initial values for x and y
    
        #Solve the system of ODEs using odeint
        solution = odeint(ODEs, initial_conditions, n, args = (lam, gam))
        main_tracks[i].set_data(solution[:,0], solution[:,1])

    #Show plots
    fig.canvas.draw()


#_________________________Cursor interactivity_______________________
#Initial cursor point
clickpoint, = ax2.plot(lam_0**2,gam_0,'r*')

#Define function for interaction with plot
def regions_plot(event):
    #Find x,y of cursor
    ix, iy = event.xdata, event.ydata

    #Update lambda, gamma values from interaction
    lambda_slide.set(np.sqrt(ix))
    gamma_slide.set(iy)
    clickpoint.set_data(ix,iy)

    #Update plot for new values
    update_plot(event)
    fig2.canvas.draw()


#____________________________Defining Widgets_____________________________________
#When cursor clicks on region plot update to clicked value
cid = fig2.canvas.mpl_connect('button_press_event', regions_plot)

#Lambda Slider
lambda_slide_label = tk.Label(window, text = 'r$\lambda$ value', width = 15, 
                       height = 2).pack(anchor = 'center', pady=0.1)
lambda_slide = tk.Scale(window, from_ = 0, to = np.sqrt(12), orient = 'horizontal',
                       width = 20, length = 250, resolution=0.01)
lambda_slide.bind("<ButtonRelease-1>", update_plot)
lambda_slide.set(lam_0)
lambda_slide.pack(anchor = 'center', pady=0.1)

#Gamma slider
gamma_slide_label = tk.Label(window, text = 'r$\gamma$ value', width = 15, 
                       height = 2).pack(anchor = 'center', pady=0.1)
gamma_slide = tk.Scale(window, from_ = 0, to = 2, orient = 'horizontal',
                       width = 20, length = 250, resolution=0.01)
gamma_slide.bind("<ButtonRelease-1>", update_plot)
gamma_slide.set(gam_0)
gamma_slide.pack(anchor = 'center', pady=0.1)


#___________________________________Initial Plot____________________________________

#Define log a linspace for ode
n = np.linspace(0, 30, 1000)
xinit = np.linspace(-0.99, 0.99, pathnum)

#Initial plot
main_tracks = []
for i in range(pathnum):
    initial_conditions = [xinit[i], 0.01]  #Initial values for x and y
    
    lam = lambda_slide.get()
    gam = gamma_slide.get()

    # Solve the system of ODEs using odeint
    solution = odeint(ODEs, initial_conditions, n, args = (lam, gam))

    solution_x = solution[:, 0]
    solution_y = solution[:, 1]

    solution_acceleration = acceleration(solution_x, solution_y)
    
    track_i = ax.plot(solution_x, solution_y, 'k', linewidth=.5)[0]
    
    main_tracks.append(track_i)
    track_i.set_visible(True)


window.mainloop()
#End
