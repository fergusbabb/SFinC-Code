#_______________________________Imports________________________________________
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
from mpl_toolkits.mplot3d import Axes3D

#Tkinter for widgets and interactivity. Scipy for ode solving
import tkinter as tk
from scipy.integrate import odeint

#This import stops automatically sets up the window so its reasonable
import ctypes 
ctypes.windll.shcore.SetProcessDpiAwareness(1) 


#Personal plotting preferences
plt.rcParams.update({"text.usetex": True, "font.family": "serif",
                     "font.serif": ["Computer Modern Serif"]})
plt.rc('axes', labelsize=12, titlesize=15)
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

#_________________________Set up main window___________________________________
#Standard tkinter window set up
window = tk.Tk()
window.title('Autonomous Systems')
window.geometry('1600x1000')

#Tracks plot figure
fig = Figure(figsize=(5, 3))
axis_dims = [.15,.20,.7,.7]
ax = fig.add_axes(axis_dims, projection='3d')
ax.set_aspect('equal')

'''
#Interactivity phase plot figure
fig2 = Figure(figsize=(5, 3))
axis_dims2 = [.15,.20,.7,.7]
ax2 = fig2.add_axes(axis_dims2)
#ax2.set_aspect('equal')

#Acceleration vs N plot figure
fig3 = Figure(figsize=(5, 3))
axis_dims3 = [.15,.20,.7,.7]
ax3 = fig3.add_axes(axis_dims3)
#ax3.set_aspect('equal')

#Acceleration vs x/y plot figure
fig4 = Figure(figsize=(5, 3))
axis_dims4 = [.15,.20,.7,.7]
ax4 = fig4.add_axes(axis_dims4)
#ax4.set_aspect('equal')
'''

#Bounding Circle
'''
theta = np.linspace(0, np.pi, 150)
ax.plot(np.cos(theta), np.sin(theta),
         'r', linewidth=1, label='Bounding Circle')
ax.plot([-1,1], [0,0], 'r', linewidth=1)
'''

#Initial values
lam_0 = 1
pathnum = 1


#_________________________Initialise plots___________________________________
#Tkinter works as:
#1) Main Window 
#2) Set up canvases inside window 
#3) Set up Figures attached to canvases 
#4) Add axes to figures as usual

ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$z$')

#plot_label = tk.Label(window, text = "Master Plot", width = 10, 
#                       height = 2).pack(anchor='n')

'''
#Interactive plot, plot lines as in P40
ax2.set(xlim=[0,12],ylim=[0,2])
ax2.set_xlabel('$\lambda^2$', x=1.02)
ax2.set_ylabel('$\gamma$', rotation = 0, y=1.02)
ax2.plot([0,6], [12,2],'k', linestyle = '-')
ax2.plot([0,12],[0,0],'k')
ax2.plot([0,0],[0,2],'k')
ax2.plot([6,6],[0,2],'k')
ax2.plot([0,6],[0,2],'k')
ax2.plot([2,2],[0,2],'k', linestyle = ':')
ax2.plot([0, 12], [2/3, 2/3],'k', linestyle = ':')
'''


canvas = FigureCanvasTkAgg(fig, window) 
#Canvas is where figure is placed to window
canvas.draw() #Show canvas (ie show figure)

'''
#As before for interactive phase plot
canvas2 = FigureCanvasTkAgg(fig2, window)
canvas2.draw()

#As before for acceleration vs N plot
canvas3 = FigureCanvasTkAgg(fig3, window)
canvas3.draw()

#As before for x/y vs N plot
canvas4 = FigureCanvasTkAgg(fig4, window)
canvas4.draw()'''


#_____________________________Define ODE for 1 fluid_____________________
s = np.sqrt(6)
def ODEs(coords, t, lam):
    # coords[0] corresponds to the variable x
    # coords[1] corresponds to the variable y
    # coords[2] corresponds to the variable y
    x, y, z = coords[0], coords[1], coords[2]

    a = 1/2
    A = a * (3*x**2 - 3*y**2 + z**2)
    
    #Return the derivatives [x', y', z']
    dx_dN = -3*a  +  3 * lam * y**2  / (x * s)  +  A
    dy_dN = -3*a  -      lam * s * x / 2        +  A
    dz_dN = -  a  +  A 

    return [dx_dN, dy_dN, dz_dN]

#____________________________Acceleration Region___________________________
'''
def acceleration(x, y, gam):
    ## Defines the acceleration expression for a given point (x, y)
    return -3*x**2 - 1.5*gam*(1 - x**2 - y**2) + 1

def acceleration_y_boundary(x, gam):
    ## Defines the line where the acceleration expression = 0 at y for a given x
    ysquared = (2/gam - 1) * x**2 - 2/(3*gam) + 1
    ysquared = np.where(ysquared>0, ysquared, 0)
    return np.sqrt(ysquared)

def compute_fill(gam):
    ## Calculates where the acceleration region is valid and draws it on the plot
    x = np.linspace(-1, 1, 401)
    boundary = acceleration_y_boundary(x, gam)
    xWhere = np.where(boundary < np.sqrt(1-x**2), True, False)
    return ax.fill_between(x, boundary, np.sqrt(1-x**2), where = xWhere)
'''  
    
#_______________________________Update track plots___________________________
def update_plot(event):
    #Before update_plot is called lam and gam sliders are updated
    #Sliders will always have current lambda/ gamma values
    lam = lambda_slide.get()
    #gam =  gamma_slide.get()

    #Update cursor star point
    #clickpoint.set_data(lam**2, gam)

y_0 =
    #Plot all paths with updated values
    for i in range(pathnum):
        z_0 = 0.9999
        y_0 = 0.0001

        initial_conditions = [xinit[i], y_0, z_0]  
        #Initial values for x y and z

        
        #Solve the system of ODEs using odeint
        solution = odeint(ODEs, initial_conditions, n, args = (lam,))

        solution_x = solution[:, 0]
        solution_y = solution[:, 1]
        solution_z = solution[:, 2]

        #solution_acceleration = acceleration(solution_x, solution_y, gam)

        main_tracks[i].set_data(solution_x, solution_y)
        main_tracks[i].set_3d_properties(solution_z)
        #acceleration_plot_tracks[i].set_ydata(solution_acceleration)
        #x_tracks[i].set_ydata(solution_x)
        #y_tracks[i].set_ydata(solution_y)
    '''
    global fill
    fill.remove()
    fill = compute_fill(gam)# = plt.fill_between(x, acceleration_y_boundary(x), np.sqrt(1-x**2), where = xWhere)
    fill.set_color([0.1, 0.4, 0.7, 0.5])
    '''
    #Show plots
    fig.canvas.draw()
    #fig2.canvas.draw()
    #fig3.canvas.draw()
    #fig4.canvas.draw()

#_________________________Cursor interactivity_______________________
#Initial cursor point
#clickpoint, = ax2.plot(lam_0**2, gam_0, 'r*')

#Define function for interaction with plot
'''def regions_plot(event):
    #Find x,y of cursor
    ix, iy = event.xdata, event.ydata

    #Update lambda, gamma values from interaction
    lambda_slide.set(np.sqrt(ix))
    gamma_slide.set(iy)
    #clickpoint.set_data(ix, iy)

    #Update plot for new values
    update_plot(event)
'''

#____________________________Defining Widgets_____________________________________
#When cursor clicks on region plot update to clicked value
#cid = fig2.canvas.mpl_connect('button_press_event', regions_plot)

#Lambda Slider
lambda_slide_label = tk.Label(window, text = 'r$\lambda$ value', width = 15, 
                       height = 2)
lambda_slide = tk.Scale(window, from_ = 0, to = np.sqrt(12), orient = 'horizontal',
                       width = 20, length = 250, resolution=0.01)
lambda_slide.bind("<ButtonRelease-1>", update_plot)
lambda_slide.set(lam_0)


#Gamma slider
'''
gamma_slide_label = tk.Label(window, text = 'r$\gamma$ value', width = 15, 
                       height = 2)
gamma_slide = tk.Scale(window, from_ = 0, to = 2, orient = 'horizontal',
                       width = 20, length = 250, resolution=0.01)
gamma_slide.bind("<ButtonRelease-1>", update_plot)
gamma_slide.set(gam_0)
'''

canvas.get_tk_widget().grid( row=1, column=1, rowspan=3, columnspan=2, ipadx=10, ipady=10)
#canvas2.get_tk_widget().grid(row=4, column=1, rowspan=3, columnspan=2, ipadx=10, ipady=10)
#canvas3.get_tk_widget().grid(row=1, column=4, rowspan=3, columnspan=2, ipadx=10, ipady=10)
#canvas4.get_tk_widget().grid(row=4, column=4, rowspan=3, columnspan=2, ipadx=10, ipady=10)
lambda_slide_label.grid(     row=1, column=3, rowspan=1, columnspan=1                    )
lambda_slide.grid(           row=1, column=3, rowspan=1, columnspan=1                    )
#gamma_slide_label.grid(      row=2, column=3, rowspan=1, columnspan=1                    )
#gamma_slide.grid(            row=2, column=3, rowspan=1, columnspan=1                    )




#___________________________________Initial Plot____________________________________

#Define log a linspace for ode
n = np.linspace(0, 10, 1000)
xinit = np.linspace(-0.99, 0.99, pathnum)

#Initial plot
main_tracks = []
acceleration_plot_tracks = []
x_tracks = []
y_tracks = []

#fill = compute_fill(gam_0)
#fill.set_color([0.1, 0.4, 0.7, 0.5])

for i in range(pathnum):
    z_0 = 0.9999
    y_0 = 0.0001

    initial_conditions = [xinit[i], y_0, z_0]  
    #Initial values for x y and z
    
    lam = lambda_slide.get()
    #gam = gamma_slide.get()

    # Solve the system of ODEs using odeint
    solution = odeint(ODEs, initial_conditions, n, args = (lam,))

    solution_x = solution[:, 0]
    solution_y = solution[:, 1]
    solution_z = solution[:, 2]
    '''
    x_plot_i = ax4.plot(n, solution_x, 'm', linewidth=.5)[0]
    x_tracks.append(x_plot_i)
    x_plot_i.set_visible(True)

    y_plot_i = ax4.plot(n, solution_y, 'c', linewidth=.5)[0]
    y_tracks.append(y_plot_i)
    y_plot_i.set_visible(True)

    solution_acceleration = acceleration(solution_x, solution_y, gam)
    
    acceleration_plot_i = ax3.plot(n, solution_acceleration, 'r', linewidth=.5)[0]
    acceleration_plot_tracks.append(acceleration_plot_i)
    acceleration_plot_i.set_visible(True)
    '''
    track_i = ax.plot(solution_x, solution_y, solution_z, 'k', linewidth=.5)[0]
    main_tracks.append(track_i)
    track_i.set_visible(True)


window.mainloop()
#End
