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
fig = Figure(figsize=(16, 10)) #1600x1000 pixels
track_axis_dims = [.05,.65,.3,.3]
track_ax = fig.add_axes(track_axis_dims, projection='3d')
track_ax.set_box_aspect([2, 1, 1])

dens_axis_dims = [.05,.15,.3,.3]
dens_ax = fig.add_axes(dens_axis_dims)

accel_axis_dims = [.45,.65,.3,.3]
accel_ax = fig.add_axes(accel_axis_dims)

gamma_axis_dims = [.45,.15,.3,.3]
gamma_ax = fig.add_axes(gamma_axis_dims)

#Bounding Circle
theta = np.linspace(0, np.pi, 150)
track_ax.plot(np.cos(theta), np.sin(theta),
         'r', linewidth=1, label='Bounding Circle')
track_ax.plot(np.cos(theta), np.zeros(150), np.sin(theta),
         'r', linewidth=1, label='Bounding Circle')
track_ax.plot([-1,1], [0,0], 'r.', linewidth=1)

#Initial values
lam_0 = 1
pathnum = 1


#_________________________Initialise plots___________________________________
#Tkinter works as:
#1) Main Window 
#2) Set up canvases inside window 
#3) Set up Figures attached to canvases 
#4) Add axes to figures as usual

track_ax.set_xlabel('$x$')
track_ax.set_ylabel('$y$')
track_ax.set_zlabel('$z$')


#_____________________________Define ODE for 2 fluids_____________________
s = np.sqrt(6)
def ODEs(state, N, lam):
    # coords[0] corresponds to the variable x
    # coords[1] corresponds to the variable y
    # coords[2] corresponds to the variable y

    x = state[0]
    y = state[1]
    z = state[2]

    HdotSection = 1 + x**2 - y**2 + (1/3)*z**2

    xprime = (-3*x) + (lam * np.sqrt(1.5) * y**2)  + (1.5*x * HdotSection)
    yprime =         (-lam * np.sqrt(1.5) * x * y) + (1.5*y * HdotSection)
    zprime = (-2*z)                                + (1.5*z * HdotSection)

    return [xprime, yprime, zprime]

#___________________________Function for fixed points________________________
def fixedPoints_func(lam):
    fixedPoints = np.array([
    [0, 0, 0],
    [0, 0, 1],
    [1, 0, 0],
    [-1, 0, 0]
    ])

    if lam**2 < 6:
        fixedPoints = np.append(fixedPoints, [[lam/np.sqrt(6), np.sqrt(1 - (lam**2)/6), 0]], axis=0)

    if lam**2 > 3:
        fixedPoints = np.append(fixedPoints, [[np.sqrt(3/2)/lam, np.sqrt(3/2)/lam, 0]], axis=0)

    if lam**2 > 8/3:
        fixedPoints = np.append(fixedPoints, [[2 * np.sqrt(2/3) / lam, 2 / (lam * np.sqrt(3)), np.sqrt(1 - (4/lam**2))]], axis=0)
    
    return fixedPoints


#_______________________________Acceleration Expression______________________
def accelerationExpression(x, y, z):
    HdotSection = 1 + x**2 - y**2 + (1/3)*z**2
    return (1/2) * (-3*x**2 + 3*y**2 - z**2 - 1)


#_______________________________Update track plots___________________________

def update_plot(event):
    #Before update_plot is called lam and gam sliders are updated
    #Sliders will always have current lambda/ gamma values
    lam = lambda_slide.get()
    fixedPoints = fixedPoints_func(lam).transpose()

    fixedPoints_plot.set_data(fixedPoints[0,:], fixedPoints[1,:])
    fixedPoints_plot.set_3d_properties(fixedPoints[2,:]) #Plot all paths with updated values
    
    for i in range(pathnum):
        z_0 = 0.9999
        y_0 = 0.0001
        x_0 = 0.001
        initial_conditions = [x_0, y_0, z_0]  
        #Initial values for x y and z

        
        # Solve the system of ODEs using odeint
        solution = odeint(ODEs, initial_conditions, n, args = (lam,))

        solution_x = solution[:, 0]
        solution_y = solution[:, 1]
        solution_z = solution[:, 2]

        main_tracks[i].set_data(solution_x, solution_y)
        main_tracks[i].set_3d_properties(solution_z)

        accel_plot.set_ydata(accelerationExpression(solution_x,solution_y,solution_z))

        Radn_dens_plot.set_ydata(solution_z**2)
        Mass_dens_plot.set_ydata(
            1 - solution_x**2 - solution_y**2 - solution_z**2)
        Phi_dens_plot.set_ydata(solution_x**2 + solution_y**2)

    #Show plots
    fig.canvas.draw()


#____________________________Defining Widgets_____________________________________
#When cursor clicks on region plot update to clicked value
#cid = fig2.canvas.mpl_connect('button_press_event', regions_plot)

canvas = FigureCanvasTkAgg(fig, window) 
#Canvas is where figure is placed to window
canvas.draw() #Show canvas (ie show figure)

#Lambda Slider
lambda_slide_label = tk.Label(window, text = '$\lambda$ value', width = 15, 
                       height = 2)
lambda_slide = tk.Scale(window, from_ = 0, to = np.sqrt(12), orient = 'horizontal',
                       width = 20, length = 250, resolution=0.01)
lambda_slide.bind("<ButtonRelease-1>", update_plot)
lambda_slide.set(lam_0)

canvas.get_tk_widget().place(relheight=1,relwidth=1)
lambda_slide_label.place(relx=0.05, rely=0.4, relheight=0.025, relwidth=0.3)
lambda_slide.place(relx=0.05, rely=0.425, relheight=0.05, relwidth=0.3)




#___________________________________Initial Plot____________________________________

#Define log a linspace for ode
n = np.linspace(0, 10, 1000)
xinit = np.linspace(-0.99, 0.99, pathnum)

#Initial plot
main_tracks = []
fixedPoints_plot, = track_ax.plot([], [], [], 'rx')

for i in range(pathnum):
    z_0 = 0.9999
    y_0 = 0.0001
    x_0 = 0.0001

    initial_conditions = [x_0, y_0, z_0]  
    #Initial values for x y and z
    
    lam = lambda_slide.get()
    fixedPoints = fixedPoints_func(lam).transpose()
    fixedPoints_plot.set_data(fixedPoints[0,:], fixedPoints[1,:])
    fixedPoints_plot.set_3d_properties(fixedPoints[2,:])

    # Solve the system of ODEs using odeint
    solution = odeint(ODEs, initial_conditions, n, args = (lam,))

    solution_x = solution[:, 0]
    solution_y = solution[:, 1]
    solution_z = solution[:, 2]

    Radn_dens_plot, = dens_ax.plot(n, solution_z**2, 'r',
            label = "$\Omega_r = z^2$")
    Mass_dens_plot, = dens_ax.plot(n,
         1 - solution_x**2 - solution_y**2 - solution_z**2, 'g',
            label = "$\Omega_m = 1 - x^2 - y^2 - z^2$")
    Phi_dens_plot, = dens_ax.plot(n, solution_x**2 + solution_y**2, 'b',
            label = "$\Omega_{\phi} = x^2 + y^2$")


    track_ax.plot(x_0,y_0,z_0, 'cx')
    track_i = track_ax.plot(solution_x, solution_y, solution_z, 'm', linewidth=2)[0]
    accel_plot, = accel_ax.plot(n, accelerationExpression(solution_x,solution_y,solution_z))

    main_tracks.append(track_i)
    track_i.set_visible(True)


dens_ax.legend()

window.mainloop()
#End
