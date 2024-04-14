#_______________________________Imports______________________________________
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
from mpl_toolkits.mplot3d import Axes3D

#Tkinter for widgets and interactivity. Scipy for ode solving
import tkinter as tk
from scipy.integrate import (odeint, quad)


#This import stops automatically sets up the window so its reasonable
import ctypes 
ctypes.windll.shcore.SetProcessDpiAwareness(1) 

#Personal plotting preferences
plt.style.use('classic')
plt.rcParams.update({"text.usetex": True, "font.family": "serif",
                     "font.serif": ["Computer Modern Serif"]})
plt.rc('axes', labelsize=12, titlesize=15)
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10



#_________________________Set up main window________________________________
#Standard tkinter window set up
window = tk.Tk()
window.title('Autonomous Systems')
window.geometry('1600x950')


#Tracks plot figure
fig = Figure(figsize=(16, 9.5)) #1600x950 pixels
fig.set_facecolor('white')
track_axis_dims = [.0,.5,.45,.5]
track_ax = fig.add_axes(track_axis_dims, projection='3d')
track_ax.set_box_aspect([2, 1, 1])
track_ax.view_init(elev=24, azim=66)

dens_axis_dims = [.55,.125,.4,.275]
dens_ax = fig.add_axes(dens_axis_dims)

accel_axis_dims = [.55,.4,.4,.275]
accel_ax = fig.add_axes(accel_axis_dims)

gamma_axis_dims = [.55,.675,.4,.275]
gamma_ax = fig.add_axes(gamma_axis_dims)

d_lum_ax_dims = [.05,.125,.35,.35] 
d_lum_ax = fig.add_axes(d_lum_ax_dims)

#Bounding Circle
theta = np.linspace(0, np.pi, 150)
track_ax.plot(np.cos(theta), np.sin(theta),
         'r', linewidth=1, label='Bounding Circle')
track_ax.plot(np.cos(theta), np.zeros(150), np.sin(theta),
         'r', linewidth=1)
track_ax.plot(np.zeros(75),np.cos(theta[0:75]),np.sin(theta[0:75]),
         'r', linewidth=1)
track_ax.plot([-1,0,1], [0,0,0], [0,0,0], 'r', linewidth=1)
track_ax.plot([0,0,0], [1,0,0], [0,0,1], 'r', linewidth=1)
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


#_____________________________Define ODE for 2 fluids_______________________
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
        fixedPoints = np.append(fixedPoints,
                                [[lam/np.sqrt(6),
                                    np.sqrt(1 - (lam**2)/6), 0]], axis=0)

    if lam**2 > 3:
        fixedPoints = np.append(fixedPoints,
                                [[np.sqrt(3/2)/lam,
                                    np.sqrt(3/2)/lam, 0]], axis=0)

    if lam**2 > 8/3:
        fixedPoints = np.append(fixedPoints,
                                [[2 * np.sqrt(2/3) / lam,
                                2 / (lam * np.sqrt(3)),
                                np.sqrt(1 - (4/lam**2))]], axis=0)
    
    return fixedPoints

#_______________________________Acceleration Expression______________________
def accelerationExpression(x, y, z):
    HdotSection = 1 + x**2 - y**2 + (1/3)*z**2
    return (1/2) * (-3*x**2 + 3*y**2 - z**2 - 1)

#_______________________________Effective Eos Parameter______________________
def gamma_phi(x, y):
    return (2*x**2) / (x**2 + y**2)

#_______________________________Integrand Form_______________________________

def integrand(N, x, y, z, N0, x0, y0, z0):
    gamma  = gamma_phi(x, y)
    O_m0   = 1 - x0**2 - y0**2 - z0**2
    O_r0   = z0**2
    O_phi0 = 1 - O_m0 - O_r0
    oneplusz = np.exp(N0 - N)

    radn_part = O_r0 * oneplusz**4
    matt_part = O_m0 * oneplusz**3
    phi_part  = O_phi0 * oneplusz**(3*gamma)
    integrand = -oneplusz/np.sqrt(radn_part + matt_part + phi_part)
    return integrand


#_______________________________Update track plots___________________________
def update_plot(event):
    #Before update_plot is called lam and gam sliders are updated
    #Sliders will always have current lambda/ gamma values
    lam = lambda_slide.get()
    fixedPoints = fixedPoints_func(lam).transpose()

    fixedPoints_plot.set_data(fixedPoints[0,:], fixedPoints[1,:])
    fixedPoints_plot.set_3d_properties(fixedPoints[2,:]) 
    #Plot all paths with updated values
    
    for i in range(pathnum):
        z_0 = 0.9999
        y_0 = 0.0001
        x_0 = 0.001
        initial_conditions = [x_0, y_0, z_0]  
        #Initial values for x y and z

        
        #Solve the system of ODEs using odeint
        solution = odeint(ODEs, initial_conditions, N, args = (lam,))
        solution_x = solution[:, 0]
        solution_y = solution[:, 1]
        solution_z = solution[:, 2]

        main_tracks[i].set_data(solution_x, solution_y)
        main_tracks[i].set_3d_properties(solution_z)

        accel_plot.set_ydata(accelerationExpression(
                    solution_x,solution_y,solution_z))

        Radn_dens_plot.set_ydata(solution_z**2)
        Mass_dens_plot.set_ydata(
            1 - solution_x**2 - solution_y**2 - solution_z**2)
        Phi_dens_plot.set_ydata(solution_x**2 + solution_y**2)

        effective_eos.set_ydata(gamma_phi(solution_x,solution_y))
        
        integrand_plot.set_ydata(integrand(N, solution_x, solution_y, solution_z, N0, x_0, y_0, z_0))
        #integrand_plotx.set_ydata(solution_x)
        #integrand_ploty.set_ydata(solution_y)
        #integrand_plotz.set_ydata(solution_z)

    #Show plots
    fig.canvas.draw()



#____________________________Defining Widgets________________________________
#When cursor clicks on region plot update to clicked value

canvas = FigureCanvasTkAgg(fig, window) 
NavigationToolbar2Tk(canvas, window)
#Show the navigation toolbar
#Canvas is where figure is placed to window
canvas.draw() #Show canvas (ie show figure)

#Lambda Slider
lambda_slide_label = tk.Label(window, text = '$\lambda$', width = 15, 
                       height = 2)
lambda_slide = tk.Scale(window, from_ = np.sqrt(6), to = 0,
                       width = 20, length = 250, resolution=0.1)

lambda_slide.set(lam_0)
lambda_slide.bind("<ButtonRelease-1>", update_plot)


canvas.get_tk_widget().place(relheight=1,relwidth=1)
lambda_slide_label.place(relx=0.45, rely=0.025,
                         relheight=0.025, relwidth=0.05)
lambda_slide.place(relx=0.45, rely=0.05, relheight=0.4, relwidth=0.05)
lambda_slide.configure(bg = 'white', borderwidth=0)
lambda_slide_label.configure(bg = 'white', borderwidth=0)


#___________________________________Initial Plot_____________________________
#Define log a linspace for ode
N0 = 0
N = np.linspace(N0, 15, 512)
a = np.exp(N)
zplus1 = 1/a

d_lum_ax.plot()

xinit = np.linspace(-0.99, 0.99, pathnum)

#Initial plot
main_tracks = []
fixedPoints_plot, = track_ax.plot([], [], [], 'o')

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
    solution = odeint(ODEs, initial_conditions, N, args = (lam,))
    solution_x = solution[:, 0]
    solution_y = solution[:, 1]
    solution_z = solution[:, 2]

    
    Radn_dens_plot, = dens_ax.plot(N, solution_z**2, 'r',
            label = "$\Omega_r = z^2$")
    Mass_dens_plot, = dens_ax.plot(N,
         1 - solution_x**2 - solution_y**2 - solution_z**2, 'g',
            label = "$\Omega_m = 1 - x^2 - y^2 - z^2$")
    Phi_dens_plot, = dens_ax.plot(N, solution_x**2 + solution_y**2, 'b',
            label = "$\Omega_\phi = x^2 + y^2$")


    track_ax.plot(x_0,y_0,z_0, 'cx')
    track_i = track_ax.plot(
                    solution_x, solution_y, solution_z, 'm', linewidth=2)[0]
    accel_plot, = accel_ax.plot(N,
                    accelerationExpression(solution_x,solution_y,solution_z))
    effective_eos, = gamma_ax.plot(N, gamma_phi(solution_x, solution_y), 'k',
            label = r'$\gamma_\phi = {2x^2}/{(x^2+y^2)}$')

    main_tracks.append(track_i)
    track_i.set_visible(True)


    integrand_plot, = d_lum_ax.plot(N0-N, integrand(N, solution_x, solution_y, solution_z, N0, x_0, y_0, z_0))
    #integrand_plotx, = d_lum_ax.plot(N,solution_x)
    #integrand_ploty, = d_lum_ax.plot(N,solution_y)
    #integrand_plotz, = d_lum_ax.plot(N,solution_z)


track_ax.set(xlabel='$x$', ylabel='$y$', zlabel='$z$',
             xticks = [-1, -0.5, 0, 0.5, 1],
             yticks = [0, 0.5, 1],
             zticks = [0, 0.5, 1])
track_ax.set_box_aspect([2, 1, 1])


accel_ax.set_ylabel("Acceleration")
accel_ax.tick_params(axis='x', which='both', labelbottom=False) 

gamma_ax.set_ylabel("$\gamma_\phi$")
gamma_ax.tick_params(axis='x', which='both', labelbottom=False) 
gamma_ax.legend()

dens_ax.set_ylabel("Density Parameters")
dens_ax.set_xlabel("$N$")
dens_ax.legend()

window.mainloop()
#End
