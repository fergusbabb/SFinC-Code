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
#plt.rcParams.update({"text.usetex": True, "font.family": "serif",
#                     "font.serif": ["Computer Modern Serif"]})
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

dens_axis_dims = [.55,.1,.4,.275]
dens_ax = fig.add_axes(dens_axis_dims)

accel_axis_dims = [.55,.4,.4,.275]
accel_ax = fig.add_axes(accel_axis_dims)

gamma_axis_dims = [.55,.7,.4,.275]
gamma_ax = fig.add_axes(gamma_axis_dims)

d_lum_ax_dims = [.075,.125,.35,.35] 
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
pathnum = 1

lam_0 = 0.38583349
Ni = -10
N = np.linspace(0, -20, 512)
z = np.exp(-N) - 1
c = 3e5     # Given in km/s
V = z * c   # ""
h = 0.738
H_0 = 100*h # km/(s Mpc)
xinit = np.linspace(-0.99, 0.99, pathnum)

lambdaConsistent = np.array([])

Omega_phi_0 = 0.73
gamma_phi_value = 0.013

Omega_r0 = 4.984e-5
peakOmega_rConsistent = np.array([])

x0Squared = gamma_phi_value / (2 * Omega_phi_0)
y0Squared = Omega_phi_0 - x0Squared

Omega_m0 = 1 - x0Squared - y0Squared - Omega_r0

state_0 = [np.sqrt(x0Squared),
           np.sqrt(y0Squared),
           np.sqrt(Omega_r0)]

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

#_______________________________Integrand Forms______________________________

def d_L_IntegrandConst(currentTotal, z, zaxis,
                        Omega_m0, Omega_r0, Omega_phi_0, path_gamma_phi):
    #Return 1/W(z) for constant LCDM model
    return 1/np.sqrt((1 - Omega_Lambda)*(1+z)**3 + 
                     (Omega_Lambda))

def d_L_IntegrandScalar(currentTotal, z, zaxis,
                         Omega_m0, Omega_r0, Omega_phi_0, path_gamma_phi):
    #Return 1/W(z) with variable scalar field
    index = np.abs(zaxis-z).argmin()
    return 1/np.sqrt(Omega_m0*(1+z)**3 +
                     Omega_r0*(1+z)**4 + 
                     Omega_phi_0*(1+z)**(3 * path_gamma_phi[index]))


#_______________________________Update track plots___________________________
def update_plot(event):
    #Before update_plot is called lam sliders are updated
    #Sliders will always have current lambda/ gamma values
    lam = lambda_slide.get()
    fixedPoints = fixedPoints_func(lam).transpose()

    #Update fixed points
    fixedPoints_plot.set_data(fixedPoints[0,:], fixedPoints[1,:])
    fixedPoints_plot.set_3d_properties(fixedPoints[2,:]) 
    
    #Plot all paths with updated values
    for i in range(pathnum):
        #Solve the system of ODEs using odeint
        solution = odeint(ODEs, state_0, N, args = (lam,))
        pathx = solution[:, 0]
        pathy = solution[:, 1]
        pathz = solution[:, 2]

        #Update tracks
        main_tracks[i].set_data(pathx, pathy)
        main_tracks[i].set_3d_properties(pathz)

        #Update acceleration plot
        accel_plot_new_data = accelerationExpression(
                                    pathx,pathy,pathz)
        accel_plot.set_ydata(accel_plot_new_data)

        #Update relative density plots
        Radn_dens_plot.set_ydata(pathz**2)
        Mass_dens_plot.set_ydata(
            1 - pathx**2 - pathy**2 - pathz**2)
        Phi_dens_plot.set_ydata(pathx**2 + pathy**2)
        
        #Update EoS plot
        path_gamma_phi = gamma_phi(pathx, pathy) 
        effective_eos.set_ydata(path_gamma_phi)
        
        #Update redshift plot
        d_L = (c/H_0) * (1 + z) * odeint(
        d_L_IntegrandScalar, 0, z, args=(
            z, Omega_m0, Omega_r0, Omega_phi_0, path_gamma_phi
            )).transpose()[0]
        integral_plot.set_ydata(d_L)

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
lambda_slide = tk.Scale(window, from_ = 1, to = 0,
                       width = 20, length = 250, resolution=0.0000000001)

lambda_slide.set(lam_0)
lambda_slide.bind("<ButtonRelease-1>", update_plot)


canvas.get_tk_widget().place(relheight=1,relwidth=1)
lambda_slide_label.place(relx=0.45, rely=0.025,
                         relheight=0.025, relwidth=0.05)
lambda_slide.place(relx=0.45, rely=0.05, relheight=0.4, relwidth=0.05)
lambda_slide.configure(bg = 'white', borderwidth=0)
lambda_slide_label.configure(bg = 'white', borderwidth=0)


#___________________________________Initial Plot_____________________________

#Initial plot
main_tracks = []
fixedPoints_plot, = track_ax.plot([], [], [], 'o')
for i in range(pathnum):
    lam = lambda_slide.get()
    fixedPoints = fixedPoints_func(lam).transpose()
    fixedPoints_plot.set_data(fixedPoints[0,:], fixedPoints[1,:])
    fixedPoints_plot.set_3d_properties(fixedPoints[2,:])
    print(lam)
    # Solve the system of ODEs using odeint
    solution = odeint(ODEs, state_0, N, args = (lam,))
    pathx = solution[:, 0]
    pathy = solution[:, 1]
    pathz = solution[:, 2]

    path_gamma_phi = gamma_phi(pathx, pathy)

    Radn_dens_plot, = dens_ax.plot(N, pathz**2, 'r',
            label = "$\Omega_r = z^2$")
    Mass_dens_plot, = dens_ax.plot(N,
         1 - pathx**2 - pathy**2 - pathz**2, 'g',
            label = "$\Omega_m = 1 - x^2 - y^2 - z^2$")
    Phi_dens_plot, = dens_ax.plot(N, pathx**2 + pathy**2, 'b',
            label = "$\Omega_\phi = x^2 + y^2$")

    x_i, y_i, z_i = state_0[0], state_0[1], state_0[2]

    track_ax.plot(x_i,y_i,z_i, 'cx')
    track_i = track_ax.plot(
                    pathx, pathy, pathz, 'm', linewidth=2)[0]
    accel_plot, = accel_ax.plot(N,
                    accelerationExpression(pathx,pathy,pathz))
    effective_eos, = gamma_ax.plot(N, gamma_phi(pathx, pathy), 'k',
            label = r'$\gamma_\phi = {2x^2}/{(x^2+y^2)}$')

    main_tracks.append(track_i)
    track_i.set_visible(True)

    
    d_L = (c/H_0) * (1 + z) * odeint(
        d_L_IntegrandScalar, 0, z, args=(
            z, Omega_m0, Omega_r0, Omega_phi_0, path_gamma_phi
            )).transpose()[0]
    

    integral_plot, = d_lum_ax.plot(d_L, V,
                        label = "$\Omega_{\phi 0} = $"+ str(Omega_phi_0))



#Plot the redshift plots for LCDM with different values of lambda
for Omega_Lambda in [0, 0.3, 0.7]:
    d_L = (c/H_0) * (1 + z) * odeint(
        d_L_IntegrandConst, 0, z, args=(
            z, Omega_m0, Omega_r0, Omega_Lambda, path_gamma_phi
            )).transpose()[0]
    
    d_lum_ax.plot(d_L, V, label = "$\Omega_\Lambda = $" + str(Omega_Lambda))



#_____________________Setting plot labels etc________________________________

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
#dens_ax.legend()

#d_lum_ax.plot(H_0 * d_L, d_L, "--", label = "$H_0d$")

d_lum_ax.set_ylabel("$d_L$ [Mpc]")
d_lum_ax.set_xlabel("$1+z$ [km/s]")
d_lum_ax.legend()


window.mainloop()
#End
