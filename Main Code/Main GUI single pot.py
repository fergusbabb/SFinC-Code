#_______________________________Imports______________________________________
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
window.title('GUI for Matter and Radiation')
window.geometry('1600x950')


#Tracks plot figure
fig = Figure(figsize=(16, 9.5)) #1600x950 pixels
fig.set_facecolor('white')
track_axis_dims = [.0,.525,.45,.5]
track_ax = fig.add_axes(track_axis_dims, projection='3d')
track_ax.view_init(elev=24, azim=66)

#Colour bar axes
cbar_ax_dims = [.4,.6,.015,.35]
cbar_ax = fig.add_axes(cbar_ax_dims)

#Relative Density axes
dens_axis_dims = [.55,.1,.4,.275]
dens_ax = fig.add_axes(dens_axis_dims)

#Acceleration axes
accel_axis_dims = [.55,.4,.4,.275]
accel_ax = fig.add_axes(accel_axis_dims)

#EoS Axes
gamma_axis_dims = [.55,.7,.4,.275]
gamma_ax = fig.add_axes(gamma_axis_dims)

#Hubble plot Axes
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


#__________________________Initial values_____________________________
pathnum = 1

lam_0 = 0.38583349
lam_min = 0
lam_max = 4
Ni = -20
NiForward = 10

N = np.linspace(0, Ni, 512)
NForward = np.linspace(0, NiForward, 256)

z = np.exp(-N) - 1
c = 1 #3e5     # Given in km/s
V = z * c   # ""
h = 0.738
H_0 = 100*h # km/(s Mpc)
xinit = np.linspace(-0.99, 0.99, pathnum)

lambdaConsistent = np.array([])

#Observed Values
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

#_________________________Initialise plots________________________________
#Tkinter works as:
#1) Main Window 
#2) Set up canvases inside window 
#3) Set up Figures attached to canvases 
#4) Add axes to figures as usual

track_ax.set_xlabel('$x$')
track_ax.set_ylabel('$y$')
track_ax.set_zlabel('$z$')


#_____________________________Define ODE for 2 fluids_______________________

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
    fixedPoints_labels = ['$O$', '$B$', '$A^+$', '$A^-$']

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
        fixedPoints_labels.append('$C$')

    if lam**2 > 3:
        fixedPoints = np.append(fixedPoints,
                                [[np.sqrt(3/2)/lam,
                                    np.sqrt(3/2)/lam, 0]], axis=0)
        fixedPoints_labels.append('$D$')

    if lam**2 > 8/3:
        fixedPoints = np.append(fixedPoints,
                                [[2 * np.sqrt(2/3) / lam,
                                2 / (lam * np.sqrt(3)),
                                np.sqrt(1 - (4/lam**2))]], axis=0)
        fixedPoints_labels.append('$E$')
    return fixedPoints, fixedPoints_labels

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
    state_0 = [float(x_entry_val.get()), 
               float(y_entry_val.get()), 
               float(z_entry_val.get())]

    for plot in fixedPoint_plots:
        plot.remove()
    fixedPoint_plots.clear()  # Clear the list of plot objects

    # Get new fixed points and plot them
    fixedPoints, fixedPoints_labels = fixedPoints_func(lam)
    for i, point in enumerate(fixedPoints):
        plot, = track_ax.plot(point[0], point[1], point[2], 'or')
        fixedPoint_plots.append(plot)

    x_i, y_i, z_i = state_0[0], state_0[1], state_0[2]
    state0_point.set_data(x_i, y_i)
    state0_point.set_3d_properties(z_i)
    
    #Plot all paths with updated values
    for i in range(pathnum):

        #Update the quiver vectors
        quiver_vectors = np.array([ODEs([pt[0], pt[1], pt[2]], N, lam)
                                    for pt in filtered_pts])
        u, v, w = quiver_vectors[:, 0], quiver_vectors[:, 1], quiver_vectors[:, 2]
        
        #3D Quiver doesn't have a built in update directions, so replot
        global quiver
        quiver.remove()
        quiver = track_ax.quiver(x_ins, y_ins, z_ins, u, v, w,
                    normalize=True, cmap=cmap, length = 0.075,
                    color=cmap(norm(magnitude)), norm=norm,
                    alpha = 0.5)
        
        #Solve the system of ODEs using odeint
        solution = odeint(ODEs, state_0, N, args = (lam,))
        pathx = solution[:, 0]
        pathy = solution[:, 1]
        pathz = solution[:, 2]
        
        #Add both ways
        solutionForward = odeint(ODEs, state_0, NForward, args = (lam,))
        pathx = np.append(pathx[::-1], solutionForward[:, 0])
        pathy = np.append(pathy[::-1], solutionForward[:, 1])
        pathz = np.append(pathz[::-1], solutionForward[:, 2])

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
            zAxis, Omega_m0, Omega_r0, Omega_phi_0, path_gamma_phi
            )).transpose()[0]
        
        integral_plot.set_xdata(d_L)

    #Show plots
    fig.canvas.draw()



#____________________________Defining Widgets________________________________
#When cursor clicks on region plot update to clicked value

#Canvas is where figure is placed to window
canvas = FigureCanvasTkAgg(fig, window)
canvas.draw() #Show canvas (ie show figure)

#Show the navigation toolbar
NavigationToolbar2Tk(canvas, window)


#Lambda Slider Label. Initialise, place, and hide weird borders
lambda_label_ax = fig.add_axes([.455,.9525,.05,.05])
lambda_label_ax.text(0,0,'$\lambda$ value')
lambda_label_ax.set_axis_off()

#Lambda Slider. Initialise, add interaction, place, hide borders
lambda_slide = tk.Scale(window, from_ = lam_min, to = lam_max,
                       width = 20, length = 250, resolution=0.001)
lambda_slide.set(lam_0)
lambda_slide.bind("<ButtonRelease-1>", update_plot)
lambda_slide.place(relx=0.45, rely=0.05, relheight=0.35, relwidth=0.05)
lambda_slide.configure(bg = 'white', borderwidth=0)


def submit():
    x0 = float(x_entry_val.get())
    y0 = float(y_entry_val.get())
    z0 = float(z_entry_val.get())

    if x0**2 + y0**2 + z0**2 <= 1:
        update_plot(0)
    else:
        print('Not valid values, check $x^2+y^2+z^2 <= 1$')

x_entry_val=tk.StringVar()
y_entry_val=tk.StringVar()
z_entry_val=tk.StringVar()

x_entry_label = tk.Label(window, text = '$x$:')
y_entry_label = tk.Label(window, text = '$y$:')
z_entry_label = tk.Label(window, text = '$z$:')

x_entry = tk.Entry(window, textvariable = x_entry_val
            ).place(relx=0.05, rely=0.45, relheight=0.05, relwidth=0.075)
y_entry = tk.Entry(window, textvariable = y_entry_val
            ).place(relx=0.15, rely=0.45, relheight=0.05, relwidth=0.075)
z_entry = tk.Entry(window, textvariable = z_entry_val
            ).place(relx=0.25, rely=0.45, relheight=0.05, relwidth=0.075)

x_entry_val.set(state_0[0])    
y_entry_val.set(state_0[1])
z_entry_val.set(state_0[2])

sub_btn=tk.Button(window, text = 'Submit', command = submit
            ).place(relx=0.45, rely=0.45, relheight=0.05, relwidth=0.05)

#Place Canvas
canvas.get_tk_widget().place(relheight=1,relwidth=1)

#________________________________Initialise the Quiver______________________________
#The method for placing anchor points for the Quivers was aided with ChatGPT

def fibonacci_sphere(max_radius, num_shells, base_points):
    all_points = []

    #Creating shells from the innermost to the outermost
    for shell in range(1, num_shells + 1):
        radius = max_radius*(shell/num_shells)
        #Decrease points with decreasing radius
        samples = int(base_points * (radius / max_radius)**2)
        phi = np.pi * (3. - np.sqrt(5.))  #Golden angle in rads

        for i in range(samples):
            y = 1 - (i / float(samples - 1)) * 2  #y goes from 1 to -1
            r = np.sqrt(1 - y * y)  #radius at y

            theta = phi*i  #golden angle increment

            x = np.cos(theta) * r * radius
            z = np.sin(theta) * r * radius
            y *= radius

            all_points.append([x, y, z])

    return np.array(all_points)

#Parameters
max_radius = 1
num_shells = 10  #Number of concentric shells
base_points = 1000  #Points on the outermost shell

#Aid of ChatGPT ends.

#Generate points nicely spaced based on fibbonaci. Makes fewer points 
#closer to the center so it looks nicely uniform 
all_points = fibonacci_sphere(max_radius, num_shells, base_points)

#Only take the points we want, add the origin
filtered_pts = np.vstack((
    all_points[
    (all_points[:, 0] >= -1) & (all_points[:, 0] <= 1) &
    (all_points[:, 1] >=  0) & (all_points[:, 1] <= 1) &
    (all_points[:, 2] >=  0) & (all_points[:, 2] <= 1)
],
[0,0,0]
))

#Finds gradient at all anchor points, sets them as the direction vectors
quiver_vectors = np.array([ODEs([pt[0], pt[1], pt[2]], N, lam_0)
                            for pt in filtered_pts])
u, v, w = quiver_vectors[:, 0], quiver_vectors[:, 1], quiver_vectors[:, 2]
x_ins, y_ins, z_ins = filtered_pts[:, 0], filtered_pts[:, 1], filtered_pts[:, 2]

#Scale for the colourbar
magnitude = np.sqrt(u**2 + v**2 + w**2).flatten()
norm = Normalize()
norm.autoscale(magnitude)
cmap = cm.spring

#Plot the Quivers
quiver = track_ax.quiver(x_ins, y_ins, z_ins, u, v, w,
                    normalize=True, cmap=cmap, length = 0.075,
                    color=cmap(norm(magnitude)), norm=norm,
                    alpha = 0.75)

#Plot the colourbar
cbar = plt.colorbar(quiver, cax=cbar_ax, orientation='vertical')
cbar.set_label('Magnitude of derivatives')


#___________________________________Initial Plot_____________________________

#Initial plot
main_tracks = []
fixedPoint_plots = []

gamma_ax.plot([N[-1], NForward[-1]], [4/3, 4/3], "k--", linewidth = 0.5)
gamma_ax.plot([N[-1], NForward[-1]], [1, 1], "k--", linewidth = 0.5)

for i in range(pathnum):
    lam = lambda_slide.get()

    fixedPoints, fixedPoints_labels = fixedPoints_func(lam_0)
    for point in fixedPoints:
        plot, = track_ax.plot(point[0], point[1], point[2], 'or')
        fixedPoint_plots.append(plot)

    # Solve the system of ODEs using odeint
    solution = odeint(ODEs, state_0, N, args = (lam,))
    pathx = solution[:, 0]
    pathy = solution[:, 1]
    pathz = solution[:, 2]

    solutionForward = odeint(ODEs, state_0, NForward, args = (lam,))
    pathx = np.append(pathx[::-1], solutionForward[:, 0])
    pathy = np.append(pathy[::-1], solutionForward[:, 1])
    pathz = np.append(pathz[::-1], solutionForward[:, 2])
    NAxis = np.append(N[::-1], NForward)
    zAxis = np.append(z[::-1], NForward)

    path_gamma_phi = gamma_phi(pathx, pathy)

    dens_ax.plot([0,0], [-0.2,1.2], 'k--')
    Radn_dens_plot, = dens_ax.plot(NAxis, pathz**2, 'r',
            label = "$\Omega_r = z^2$")
    Mass_dens_plot, = dens_ax.plot(NAxis,
         1 - pathx**2 - pathy**2 - pathz**2, 'g',
            label = "$\Omega_m = 1 - x^2 - y^2 - z^2$")
    Phi_dens_plot, = dens_ax.plot(NAxis, pathx**2 + pathy**2, 'b',
            label = "$\Omega_\phi = x^2 + y^2$")

    x_i, y_i, z_i = state_0[0], state_0[1], state_0[2]


    track_i = track_ax.plot(
                    pathx, pathy, pathz, 'm', linewidth=2)[0]
    state0_point, = track_ax.plot(x_i,y_i,z_i, 'cX')
    accel_plot, = accel_ax.plot(NAxis,
                    accelerationExpression(pathx,pathy,pathz))
    effective_eos, = gamma_ax.plot(NAxis, gamma_phi(pathx, pathy), 'k',
            label = r'$\gamma_\phi = {2x^2}/{(x^2+y^2)}$')

    main_tracks.append(track_i)
    track_i.set_visible(True)

    
    d_L = (c/H_0) * (1 + z) * odeint(
        d_L_IntegrandScalar, 0, z, args=(
            zAxis, Omega_m0, Omega_r0, Omega_phi_0, path_gamma_phi
            )).transpose()[0]
    

    integral_plot, = d_lum_ax.plot(d_L, V,
                        label = "$\Omega_{\phi 0} = $"+ str(Omega_phi_0))



#Plot the redshift plots for LCDM with different values of Cosm. Const.
for Omega_Lambda in [0.65, 0.7, 0.75]:
    d_L = (c/H_0) * (1 + z) * odeint(
        d_L_IntegrandConst, 0, z, args=(
            zAxis, Omega_m0, Omega_r0, Omega_Lambda, path_gamma_phi
            )).transpose()[0]
    
    d_lum_ax.plot(d_L, z, label = "$\Omega_\Lambda = $" + str(Omega_Lambda))



#_____________________Setting plot labels etc________________________________

track_ax.set(xlabel='$x$', ylabel='$y$', zlabel='$z$',
             xlim = [-1,1], ylim = [0,1], zlim = [0,1],
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
d_lum_ax.set_xlabel("$z$ [km/s]")
d_lum_ax.set_xlim(left = 0.01)
#d_lum_ax.legend(loc=4)
d_lum_ax.set_xscale('log', base=10, subs=[10**x
                         for x in (0.25, 0.5, 0.75)], nonpositive='mask')
d_lum_ax.set_yscale('log', base=10, subs=[10**x
                         for x in (0.25, 0.5, 0.75)], nonpositive='mask')



#Run the code
window.mainloop()
#End
