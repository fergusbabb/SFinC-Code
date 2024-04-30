#_______________________________Imports________________________________________
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
from matplotlib.colors import Normalize
import matplotlib.cm as cm


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
window.geometry('1600x950')


#Tracks plot figure
fig = Figure(figsize=(16, 9.5)) #1600x950 pixels
fig.set_facecolor('white')
track_axis_dims = [.0,.5,.45,.5]
track_ax = fig.add_axes(track_axis_dims)

#Colour bar axes
cbar_ax_dims = [.4,.6,.015,.35]
cbar_ax = fig.add_axes(cbar_ax_dims)

#Hubble plot Axes
lam_gam_dims = [.075,.125,.35,.35] 
lam_gam_ax = fig.add_axes(lam_gam_dims)

#Bounding Circle
theta = np.linspace(0, np.pi, 150)
track_ax.plot(np.cos(theta), np.sin(theta),
         'r', linewidth=1, label='Bounding Circle')
track_ax.plot([-1,1], [0,0], 'r', linewidth=1)


#Initial values
gam_0 = 1
lam_0 = 3
pathnum = 50
N = np.linspace(0, 8, 1000) 
xinit = np.linspace(-0.99, 0.99, pathnum)


#_________________________Initialise plots___________________________________
#Tkinter works as:
#1) Main Window 
#2) Set up canvases inside window 
#3) Set up Figures attached to canvases 
#4) Add axes to figures as usual

track_ax.set_xlabel('$x$', x=1.02)
track_ax.set_ylabel('$y$', rotation = 0, y=1.02)

#plot_label = tk.Label(window, text = "Master Plot", width = 10, 
#                       height = 2).pack(anchor='n')

#Interactive plot, plot lines as in P40
lam_gam_ax.set(xlim=[0,12],ylim=[0,2])
lam_gam_ax.set_xlabel('$\lambda^2$', x=1.02)
lam_gam_ax.set_ylabel('$\gamma$', rotation = 0, y=1.02)
lam_gam_ax.plot([0,6], [12,2],'k', linestyle = '-')
lam_gam_ax.plot([0,12],[0,0],'k')
lam_gam_ax.plot([0,0],[0,2],'k')
lam_gam_ax.plot([6,6],[0,2],'k')
lam_gam_ax.plot([0,6],[0,2],'k')
lam_gam_ax.plot([2,2],[0,2],'k', linestyle = ':')
lam_gam_ax.plot([0, 12], [2/3, 2/3],'k', linestyle = ':')


canvas = FigureCanvasTkAgg(fig, window) 
#Canvas is where figure is placed to window
canvas.draw() #Show canvas (ie show figure)

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

#____________________________Acceleration Region___________________________
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
    return track_ax.fill_between(x, boundary, np.sqrt(1-x**2), where = xWhere)
    
    
#_______________________________Update track plots___________________________
def update_plot(event):
    #Before update_plot is called lam and gam sliders are updated
    #Sliders will always have current lambda/ gamma values
    lam = lambda_slide.get()
    gam =  gamma_slide.get()

    #Update cursor star point
    clickpoint.set_data(lam**2, gam)

    #Plot all paths with updated values
    for i in range(pathnum):
        initial_conditions = [xinit[i], 0.01]  #Initial values for x and y

        
        #Solve the system of ODEs using odeint
        solution = odeint(ODEs, initial_conditions, n, args = (lam, gam))

        solution_x = solution[:, 0]
        solution_y = solution[:, 1]

        solution_acceleration = acceleration(solution_x, solution_y, gam)

        main_tracks[i].set_data(solution_x, solution_y)
        acceleration_plot_tracks[i].set_ydata(solution_acceleration)
        x_tracks[i].set_ydata(solution_x)
        y_tracks[i].set_ydata(solution_y)

        #Update the quiver vectors
        quiver_vectors = np.array([ODEs([pt[0], pt[1]], N, lam, gam)
                                    for pt in points])
        u, v = quiver_vectors[:, 0], quiver_vectors[:, 1]
        magnitude = np.sqrt(u**2 + v**2)
        u = u/magnitude
        v = v/magnitude

        global quiver
        quiver.remove()
        quiver = track_ax.quiver(x_ins, y_ins, u, v,
                    color=cmap(norm(magnitude)), norm=norm,
                    cmap=cmap)
        
    global fill
    fill.remove()
    fill = compute_fill(gam)# = plt.fill_between(x, acceleration_y_boundary(x), np.sqrt(1-x**2), where = xWhere)
    fill.set_color([0.1, 0.4, 0.7, 0.5])

    #Show plots
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
magnitude = np.sqrt(u**2 + v**2)
norm = Normalize()
u = u/magnitude
v = v/magnitude

cmap = cm.spring

#Plot the Quivers
quiver = track_ax.quiver(x_ins, y_ins, u, v,
                    color=cmap(norm(magnitude)), norm=norm,
                    cmap=cmap, scale=25, headwidth = 3)

#Plot the colourbar
cbar = plt.colorbar(quiver, cax=cbar_ax, orientation='vertical')
cbar.set_label('Magnitude of derivatives')

#_________________________Cursor interactivity_______________________
#Initial cursor point
clickpoint, = lam_gam_ax.plot(lam_0**2, gam_0, 'r*')

#Define function for interaction with plot
def regions_plot(event):
    #Find x,y of cursor
    ix, iy = event.xdata, event.ydata

    #Update lambda, gamma values from interaction
    lambda_slide.set(np.sqrt(ix))
    gamma_slide.set(iy)
    clickpoint.set_data(ix, iy)

    #Update plot for new values
    update_plot(event)


#____________________________Defining Widgets_____________________________________
#When cursor clicks on region plot update to clicked value
cid = fig.canvas.mpl_connect('button_press_event', regions_plot)

#Lambda Slider
lambda_slide_label = tk.Label(window, text = 'r$\lambda$ value', width = 15, 
                       height = 2)
lambda_slide = tk.Scale(window, from_ = 0, to = np.sqrt(12), orient = 'horizontal',
                       width = 20, length = 250, resolution=0.01)
lambda_slide.bind("<ButtonRelease-1>", update_plot)
lambda_slide.set(lam_0)


#Gamma slider
gamma_slide_label = tk.Label(window, text = 'r$\gamma$ value', width = 15, 
                       height = 2)
gamma_slide = tk.Scale(window, from_ = 0, to = 2, orient = 'horizontal',
                       width = 20, length = 250, resolution=0.01)
gamma_slide.bind("<ButtonRelease-1>", update_plot)
gamma_slide.set(gam_0)


canvas.get_tk_widget().grid( row=1, column=1, rowspan=3, columnspan=2, ipadx=10, ipady=10)
lambda_slide_label.grid(     row=1, column=3, rowspan=1, columnspan=1                    )
lambda_slide.grid(           row=1, column=3, rowspan=1, columnspan=1                    )
gamma_slide_label.grid(      row=2, column=3, rowspan=1, columnspan=1                    )
gamma_slide.grid(            row=2, column=3, rowspan=1, columnspan=1                    )




#___________________________________Initial Plot____________________________________

#Define log a linspace for ode
n = np.linspace(0, 10, 1000)
xinit = np.linspace(-0.99, 0.99, pathnum)

#Initial plot
main_tracks = []
acceleration_plot_tracks = []
x_tracks = []
y_tracks = []

fill = compute_fill(gam_0)
fill.set_color([0.1, 0.4, 0.7, 0.5])

for i in range(pathnum):
    initial_conditions = [xinit[i], 0.01]  #Initial values for x and y
    
    lam = lambda_slide.get()
    gam = gamma_slide.get()

    # Solve the system of ODEs using odeint
    solution = odeint(ODEs, initial_conditions, n, args = (lam, gam))

    solution_x = solution[:, 0]
    solution_y = solution[:, 1]

    #x_plot_i = ax4.plot(n, solution_x, 'm', linewidth=.5)[0]
    #x_tracks.append(x_plot_i)
    #x_plot_i.set_visible(True)

    #y_plot_i = ax4.plot(n, solution_y, 'c', linewidth=.5)[0]
    #y_tracks.append(y_plot_i)
    #y_plot_i.set_visible(True)

    solution_acceleration = acceleration(solution_x, solution_y, gam)
    
    #acceleration_plot_i = ax3.plot(n, solution_acceleration, 'r', linewidth=.5)[0]
    #acceleration_plot_tracks.append(acceleration_plot_i)
    #acceleration_plot_i.set_visible(True)

    track_i = track_ax.plot(solution_x, solution_y, 'k', linewidth=.5, alpha=0.4)[0]
    main_tracks.append(track_i)
    track_i.set_visible(True)


window.mainloop()
#End
