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
ax.set_box_aspect([2, 1, 1])

#Bounding Circle

theta = np.linspace(0, np.pi, 150)
ax.plot(np.cos(theta), np.sin(theta),
         'r', linewidth=1, label='Bounding Circle')
ax.plot(np.cos(theta), np.zeros(150), np.sin(theta),
         'r', linewidth=1, label='Bounding Circle')
ax.plot([-1,1], [0,0], 'r.', linewidth=1)

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

canvas = FigureCanvasTkAgg(fig, window) 
#Canvas is where figure is placed to window
canvas.draw() #Show canvas (ie show figure)


#_____________________________Define ODE for 1 fluid_____________________
s = np.sqrt(6)
def ODEs(state, N, lam):
    # coords[0] corresponds to the variable x
    # coords[1] corresponds to the variable y
    # coords[2] corresponds to the variable y

    x = state[0]
    y = state[1]
    z = state[2]

    HdotSection = 1 + x**2 - y**2 - (1/3)*z**2

    xprime = (-3*x) + (lam * np.sqrt(1.5) * y**2)  + (1.5*x * HdotSection)
    yprime =         (-lam * np.sqrt(1.5) * x * y) + (1.5*y * HdotSection)
    zprime = (-2*z)                                + (1.5*z * HdotSection)

    return [xprime, yprime, zprime]

#_______________________________Update track plots___________________________

def update_plot(event):
    #Before update_plot is called lam and gam sliders are updated
    #Sliders will always have current lambda/ gamma values
    lam = lambda_slide.get()

    #Plot all paths with updated values
    for i in range(pathnum):
        z_0 = 0.9999
        y_0 = 0.0001
        x_0 = 0.001
        initial_conditions = [x_0, y_0, z_0]  
        #Initial values for x y and z

        
        #Solve the system of ODEs using odeint
        solution = odeint(ODEs, initial_conditions, n, args = (lam,))

        solution_x = solution[:, 0]
        solution_y = solution[:, 1]
        solution_z = solution[:, 2]

        main_tracks[i].set_data(solution_x, solution_y)
        main_tracks[i].set_3d_properties(solution_z)

    #Show plots
    fig.canvas.draw()


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


canvas.get_tk_widget().grid( row=1, column=1, rowspan=3, columnspan=2, ipadx=10, ipady=10)
lambda_slide_label.grid(     row=1, column=3, rowspan=1, columnspan=1                    )
lambda_slide.grid(           row=1, column=3, rowspan=1, columnspan=1                    )




#___________________________________Initial Plot____________________________________

#Define log a linspace for ode
n = np.linspace(0, 10, 1000)
xinit = np.linspace(-0.99, 0.99, pathnum)

#Initial plot
main_tracks = []
acceleration_plot_tracks = []
x_tracks = []
y_tracks = []

for i in range(pathnum):
    z_0 = 0.9999
    y_0 = 0.0001
    x_0 = 0.0001

    initial_conditions = [x_0, y_0, z_0]  
    #Initial values for x y and z
    
    lam = lambda_slide.get()
    #gam = gamma_slide.get()

    # Solve the system of ODEs using odeint
    solution = odeint(ODEs, initial_conditions, n, args = (lam,))

    solution_x = solution[:, 0]
    solution_y = solution[:, 1]
    solution_z = solution[:, 2]

    ax.plot(x_0,y_0,z_0, 'cx')
    track_i = ax.plot(solution_x, solution_y, solution_z, 'm', linewidth=2)[0]
    main_tracks.append(track_i)
    track_i.set_visible(True)


window.mainloop()
#End
