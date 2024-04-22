#_______________________________Imports________________________________________
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
from mpl_toolkits.mplot3d import Axes3D
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
window.geometry('750x400')

#Tracks plot figure
fig = Figure(figsize=(7, 5))
axis_dims = [.0,.1,.8,.8]
ax = fig.add_axes(axis_dims, projection='3d')
ax.set_box_aspect([2, 1, 1])
ax.view_init(elev=24, azim=66)

#Bounding Circle
theta_bind = np.linspace(0, np.pi, 150)
ax.plot(np.cos(theta_bind), np.sin(theta_bind),
         'r', linewidth=1, label='Bounding Circle')
ax.plot(np.cos(theta_bind), np.zeros(150), np.sin(theta_bind),
         'r', linewidth=1)
ax.plot(np.zeros(75),np.cos(theta_bind[0:75]),np.sin(theta_bind[0:75]),
         'r', linewidth=1)
ax.plot([-1,0,1], [0,0,0], [0,0,0], 'r', linewidth=1)
ax.plot([0,0,0], [1,0,0], [0,0,1], 'r', linewidth=1)

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

        #Update the quiver vectors
        vectors = np.array([ODEs([pt[0], pt[1], pt[2]], n, lam) for pt in filtered_points])
        u, v, w = vectors[:, 0], vectors[:, 1], vectors[:, 2]
        global quiver
        quiver.remove()
        quiver = ax.quiver(x_ins, y_ins, z_ins, u, v, w,
                    normalize=True, cmap=cmap, length = 0.075,
                    color=cmap(norm(magnitude)), norm=norm,
                    alpha = 0.5)
        
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
lambda_slide_label = tk.Label(window, text = '$\lambda$', width = 15, 
                       height = 2)
lambda_slide = tk.Scale(window, from_ = 0, to = np.sqrt(12), orient = 'vertical',
                       width = 20, length = 250, resolution=0.01)
lambda_slide.bind("<ButtonRelease-1>", update_plot)
lambda_slide.set(lam_0)

canvas.get_tk_widget().place(relheight=1,relwidth=1)
lambda_slide_label.place(relx=0.825, rely=0.025,
                         relheight=0.1, relwidth=0.15)
lambda_slide.place(relx=0.85, rely=0.1, relheight=0.75, relwidth=0.1)
lambda_slide.configure(bg = 'white', borderwidth=0)
lambda_slide_label.configure(bg = 'white', borderwidth=0)



#________________________________Initialise the Quiver______________________________

def scaled_fibonacci_sphere(max_radius, num_shells, base_points):
    all_points = []
    # Creating shells from the innermost to the outermost
    for shell in range(1, num_shells + 1):
        radius = max_radius * (shell / num_shells)
        # Decrease points quadratically with decreasing radius
        samples = int(base_points * (radius / max_radius)**2)
        phi = np.pi * (3. - np.sqrt(5.))  # Golden angle in radians

        for i in range(samples):
            y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
            r = np.sqrt(1 - y * y)  # radius at y

            theta = phi * i  # golden angle increment

            x = np.cos(theta) * r * radius
            z = np.sin(theta) * r * radius
            y *= radius

            all_points.append([x, y, z])

    return np.array(all_points)

# Parameters
max_radius = 1
num_shells = 10  # Number of concentric shells
base_points = 1000  # Points on the outermost shell

# Generate points
all_points = scaled_fibonacci_sphere(max_radius, num_shells, base_points)

filtered_points = np.vstack((
    all_points[
    (all_points[:, 0] >= -1) & (all_points[:, 0] <= 1) &
    (all_points[:, 1] >=  0) & (all_points[:, 1] <= 1) &
    (all_points[:, 2] >=  0) & (all_points[:, 2] <= 1)
],
[0,0,0]
))


#___________________________________Initial Plot____________________________________

#Define log a linspace for ode
n = np.linspace(0, 10, 1000)
xinit = np.linspace(-0.99, 0.99, pathnum)

#Initial plot
main_tracks = []
acceleration_plot_tracks = []
x_tracks = []
y_tracks = []


vectors = np.array([ODEs([pt[0], pt[1], pt[2]], n, lam_0) for pt in filtered_points])
u, v, w = vectors[:, 0], vectors[:, 1], vectors[:, 2]
x_ins, y_ins, z_ins = filtered_points[:, 0], filtered_points[:, 1], filtered_points[:, 2]

magnitude = np.sqrt(u**2 + v**2 + w**2).flatten()
norm = Normalize()
norm.autoscale(magnitude)

cmap = cm.spring
quiver = ax.quiver(x_ins, y_ins, z_ins, u, v, w,
                    normalize=True, cmap=cmap, length = 0.075,
                    color=cmap(norm(magnitude)), norm=norm,
                    alpha = 0.5)

cbar = plt.colorbar(quiver, ax=ax, orientation='vertical')
cbar.set_label('Magnitude of derivatives')

for i in range(pathnum):
    z_0 = 0.9999
    y_0 = 0.0001
    x_0 = 0.0001

    initial_conditions = [x_0, y_0, z_0]  


    lam = lambda_slide.get()
    #gam = gamma_slide.get()
    

    # Solve the system of ODEs using odeint
    solution = odeint(ODEs, initial_conditions, n, args = (lam,))

    solution_x = solution[:, 0]
    solution_y = solution[:, 1]
    solution_z = solution[:, 2]

    
    track_i = ax.plot(solution_x, solution_y, solution_z, 'm', linewidth=2)[0]
    main_tracks.append(track_i)
    track_i.set_visible(True)
    ax.plot(x_0,y_0,z_0, 'co')


window.mainloop()
#End
