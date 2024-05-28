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
from scipy.integrate import (odeint, quad, cumulative_trapezoid)


#This import stops automatically sets up the window so its reasonable
import ctypes 
ctypes.windll.shcore.SetProcessDpiAwareness(1) 

#Personal plotting preferences
plt.rcParams.update({"text.usetex": True, "font.family": "serif",
                     "font.serif": ["Computer Modern Serif"]})
plt.rc('axes', labelsize=14, titlesize=15)
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12


#_________________________Set up main window________________________________
#Main GUI figure
window_gui = tk.Tk()
window_gui.title('GUI for Matter and Radiation')
window_gui.geometry('1600x950')
fig = Figure(figsize=(16, 9.5)) #1600x950 pixels
fig.set_facecolor('white')

#Define Axes
#Tracks plot figure
track_axis_dims = [-.025,.525,.45,.5]
track_ax = fig.add_axes(track_axis_dims, projection='3d')
track_ax.view_init(elev=24, azim=66)

cbar_ax_dims = [.375,.6,.015,.35]
cbar_ax = fig.add_axes(cbar_ax_dims)

#Relative Density axes
dens_axis_dims = [.625,.125,.325,.3675]
dens_ax = fig.add_axes(dens_axis_dims)

#EoS Axes
gamma_axis_dims = [.625,.6,.325,.35]
gamma_ax = fig.add_axes(gamma_axis_dims)

#Hubble plot Axes

hubble_ax_dims = [.1,.125,.325,.3675]
hubble_ax = fig.add_axes(hubble_ax_dims)


#__________________________Initial values_____________________________

#Bounding Circle
theta = np.linspace(0, np.pi, 150)
track_ax.plot(np.cos(theta), np.sin(theta),
         'k', linewidth=1, label='Bounding Circle')
track_ax.plot(np.cos(theta), np.zeros(150), np.sin(theta),
         'k', linewidth=1)
track_ax.plot(np.zeros(75),np.cos(theta[0:75]),np.sin(theta[0:75]),
         'k', linewidth=1)
track_ax.plot([-1,0,1], [0,0,0], [0,0,0], 'k', linewidth=1)
track_ax.plot([0,0,0], [1,0,0], [0,0,1], 'k', linewidth=1)

pathnum = 1
# So we get B, B, BC, BEC,  EC
lam_0s = [0, np.sqrt(1), np.sqrt(4), np.sqrt(100)]
lam_0 = lam_0s[1]
lam_min = 0
lam_max = 10

Ni = 25

N = np.linspace(0, Ni, abs(int(Ni*10000)))


c = 1 #3e5     # Given in km/s
h = 0.738
H_0 = 100*h # km/(s Mpc)
xinit = np.linspace(-0.99, 0.99, pathnum)

H0_SN = 73.04
H0_SN_Err = 1.04

H0_PL = 68.75
H0_PL_Err = 0.52

#Observed Values
Omega_phi_0 = 0.73
gamma_phi_value = 0.013
Omega_r0 = 4.984e-5

peakOmega_rConsistent = np.array([])

x0Squared = gamma_phi_value / (2 * Omega_phi_0)
y0Squared = Omega_phi_0 - x0Squared

Omega_m0 = 1 - x0Squared - y0Squared - Omega_r0

#For lambda^2 = 1,4 we use x = 2e-8, y=1e-9, z=0.99, Nplotmax = 7
state_0 = [2e-11,
           1e-11,
           0.999]
# state_0 = [2e-8,
#            1e-9,
#            0.99]

N_plot_max = 7
N_plot_min = -15

#_________________________Initialise plots________________________________
#Tkinter works as:
#1) Main Window 
#2) Set up canvases inside window 
#3) Set up Figures attached to canvases 
#4) Add axes to figures as usual


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
    fixedPoints_labels = ['$O$', '$A^+$', '$A^-$']

    fixedPoints = np.array([
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0]
    ])

    if lam**2 < 6:
        fixedPoints = np.append(fixedPoints,
                                [[lam/np.sqrt(6),
                                    np.sqrt(1 - (lam**2)/6), 0]], axis=0)
        fixedPoints_labels.append('$B$')

    if lam**2 >3:
        fixedPoints = np.append(fixedPoints,
                                [[np.sqrt(3/2)/lam, np.sqrt(3/2)/lam, 0]], axis=0)
        fixedPoints_labels.append('$C$')

    fixedPoints = np.append(fixedPoints,
                            [[0, 0, 1]], axis=0)
    fixedPoints_labels.append('$D$')

    if lam**2 > 4:
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

#_______________________________Find Redshift from e-foldings________________
def getRedshift(N):
    return np.exp(-N) - 1

#_______________________________Effective Eos Parameter______________________
def gamma_phi(x, y):
    return (2*x**2) / (x**2 + y**2)

#_______________________________Integrand Forms______________________________

def d_L_Dark_Energy(currentTotal, z, zaxis, Omega_m0, Omega_r0, Omega_Lambda0, w_Lambda):
    return 1/np.sqrt(Omega_m0*(1+z)**3 +
                     Omega_r0*(1+z)**4 + 
                     Omega_Lambda0*(1+z)**(3 * (1 + w_Lambda)))

#def d_L_IntegrandConst(currentTotal, z, zaxis,
#                        Omega_m0, Omega_r0, Omega_phi_0, path_gamma_phi):
#    #Return 1/W(z) for constant LCDM model
#    return 1/np.sqrt((1 - Omega_Lambda)*(1+z)**3 + 
#                     (Omega_Lambda))

def d_L_IntegrandScalar(currentTotal, z, zaxis,
                         Omega_m0, Omega_r0, Omega_phi_0, path_gamma_phi):
    #Return 1/W(z) with variable scalar field
    index = np.abs(zaxis-z).argmin()
    return 1/np.sqrt(Omega_m0*(1+z)**3 +
                     Omega_r0*(1+z)**4 + 
                     Omega_phi_0*(1+z)**(3 * path_gamma_phi[index]))


#____________________________Calculate H from x and y________________________

def HfromY(pathy, pathxIntegral, lam, V0=1):
    kappa = 1#?
    exponents = np.exp(- lam * kappa * pathxIntegral / 2)
    H = kappa * np.sqrt(V0) * exponents / (pathy * np.sqrt(3))
    return H


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
        plot, = track_ax.plot(point[0], point[1], point[2], 'or', markeredgewidth=0.5, markeredgecolor='black')
        fixedPoint_plots.append(plot)

    
    x_i, y_i, z_i = state_0[0], state_0[1], state_0[2]
    state0_point.set_data(x_i, y_i)
    state0_point.set_3d_properties(z_i)
    

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
                alpha = 0.75, linewidth=1)
    
    #Solve the system of ODEs using odeint
    solution = odeint(ODEs, state_0, N, args = (lam,))
    pathx = solution[:, 0]
    pathy = solution[:, 1]
    pathz = solution[:, 2]

    #Update tracks
    track.set_data(pathx, pathy)
    track.set_3d_properties(pathz)

    #Find when in N each event occured
    NAxis = N

    rad_dens = pathz**2
    mass_dens = 1 - pathx**2 - pathy**2 - pathz**2
    phi_dens = pathx**2 + pathy**2
    total_dens = rad_dens + mass_dens + phi_dens
    
    indexToday = np.argmin(np.abs(phi_dens-Omega_phi_0))
    indexMR_eq = np.argmin(np.abs(mass_dens-rad_dens)[0:indexToday])
    indexMPhi_eq = np.argmin(np.abs(mass_dens-phi_dens)[indexMR_eq:-1]) + indexMR_eq
    indexMPeak = np.argmax(mass_dens)

    NAxis -= NAxis[indexToday]
    zAxis = getRedshift(NAxis[:indexToday + 1])
    z = zAxis[::-1]
    
    MR_eqLine.set_xdata([NAxis[indexMR_eq], NAxis[indexMR_eq]])
    MPhi_eqLine.set_xdata([NAxis[indexMPhi_eq], NAxis[indexMPhi_eq]])
    MPeakLine.set_xdata([NAxis[indexMPeak], NAxis[indexMPeak]])

    #Update relative density plots
    Radn_dens_plot.set_ydata(pathz**2)
    Radn_dens_plot.set_xdata(NAxis)
    
    Mass_dens_plot.set_ydata(1 - pathx**2 - pathy**2 - pathz**2)
    Mass_dens_plot.set_xdata(NAxis)
    
    Phi_dens_plot.set_ydata(pathx**2 + pathy**2)
    Phi_dens_plot.set_xdata(NAxis)
    
    #Update acceleration plot
    accel_plot_new_data = accelerationExpression(
                                pathx,pathy,pathz)
    accel_plot.set_ydata(accel_plot_new_data)
    accel_plot.set_xdata(NAxis)

    gamma_ax.set_xlim([NAxis[0],N_plot_max])
    dens_ax.set_xlim([NAxis[0],N_plot_max])
    
    rScalingLine.set_xdata([NAxis[0],N_plot_max])
    mScalingLine.set_xdata([NAxis[0],N_plot_max])
    CCScalingLine.set_xdata([NAxis[0],N_plot_max])

    mr_eq_val = getRedshift(NAxis[indexMR_eq])
    mr_eq_text.set_text(f'$\Omega_m=\Omega_r:\; z={mr_eq_val:.3f}$')

    m_max_val = getRedshift(NAxis[indexMPeak])
    m_max_text.set_text(f'max$(\Omega_m):\; z={m_max_val:.3f}$')

    msf_eq_val = getRedshift(NAxis[indexMPhi_eq])
    msf_eq_text.set_text(f'$\Omega_m=\Omega_\phi:\; z={msf_eq_val:.3f}$')

    #Update EoS plot
    path_gamma_phi = gamma_phi(pathx, pathy) 
    effective_eos.set_ydata(path_gamma_phi)
    effective_eos.set_xdata(NAxis)


    xRunningIntegral = np.append(0, cumulative_trapezoid(pathx, NAxis))

    hubbleFromY = HfromY(pathy, xRunningIntegral, lam)
    hubbleFromY *= H_0 / hubbleFromY[indexToday]

    # Update Hubble Plot
    hubble_plot.set_ydata(hubbleFromY)
    hubble_plot.set_xdata(NAxis)

    #Show plots
    fig.canvas.draw()


#____________________________Defining Widgets________________________________
#When cursor clicks on region plot update to clicked value

#Canvas is where figure is placed to window
canvas = FigureCanvasTkAgg(fig, window_gui)
canvas.draw() #Show canvas (ie show figure)

#Show the navigation toolbar
NavigationToolbar2Tk(canvas, window_gui)


#Lambda Slider Label. Initialise, place, and hide weird borders
lambda_label_ax = fig.add_axes([.455,.9525,.05,.05])
lambda_label_ax.text(0,0,'$\lambda$ value')
lambda_label_ax.set_axis_off()

#Lambda Slider. Initialise, add interaction, place, hide borders
lambda_slide = tk.Scale(window_gui, from_ = lam_min, to = lam_max,
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

x_entry_label_ax = fig.add_axes([0.075,.525,.05,.075])
x_entry_label_ax.text(0,0,'$x_0$:')
x_entry_label_ax.set_axis_off()

y_entry_label_ax = fig.add_axes([0.175,.525,.05,.075])
y_entry_label_ax.text(0,0,'$y_0$:')
y_entry_label_ax.set_axis_off()

z_entry_label_ax = fig.add_axes([0.275,.525,.05,.075])
z_entry_label_ax.text(0,0,'$z_0$:')
z_entry_label_ax.set_axis_off()


x_entry = tk.Entry(window_gui, textvariable = x_entry_val, relief='solid'
            ).place(relx=0.1, rely=0.457, relheight=0.03, relwidth=0.05)
y_entry = tk.Entry(window_gui, textvariable = y_entry_val, relief='solid'
            ).place(relx=0.2, rely=0.457, relheight=0.03, relwidth=0.05)
z_entry = tk.Entry(window_gui, textvariable = z_entry_val, relief='solid'
            ).place(relx=0.3, rely=0.457, relheight=0.03, relwidth=0.05)

x_entry_val.set(state_0[0])    
y_entry_val.set(state_0[1])
z_entry_val.set(state_0[2])

sub_btn=tk.Button(window_gui, text = 'Submit', command = submit
            ).place(relx=0.375, rely=0.45, relheight=0.05, relwidth=0.05)

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
                    alpha = 0.75, linewidth=1)

#Plot the colourbar
cbar = plt.colorbar(quiver, cax=cbar_ax, orientation='vertical')
cbar.set_label('Magnitude of derivatives')


#___________________________________Initial Plot_____________________________

#Initial plot
fixedPoint_plots = []



lam = lambda_slide.get()

# Solve the system of ODEs using odeint
solution = odeint(ODEs, state_0, N, args = (lam,))
pathx = solution[:, 0]
pathy = solution[:, 1]
pathz = solution[:, 2]

path_gamma_phi = gamma_phi(pathx, pathy)

NAxis = N   

rad_dens = pathz**2
mass_dens = 1 - pathx**2 - pathy**2 - pathz**2
phi_dens = pathx**2 + pathy**2
total_dens = rad_dens + mass_dens + phi_dens

indexToday = np.argmin(np.abs(phi_dens-Omega_phi_0))


#Find when in N each event occured
indexMR_eq = np.argmin(np.abs(mass_dens-rad_dens)[0:indexToday])
indexMPhi_eq = np.argmin(np.abs(mass_dens-phi_dens)[indexMR_eq:-1]) + indexMR_eq
indexMPeak = np.argmax(mass_dens)

NAxis -= NAxis[indexToday]
zAxis = getRedshift(NAxis[:indexToday + 1])
z = zAxis[::-1]

mr_eq_val = getRedshift(NAxis[indexMR_eq])
mr_eq_ax = fig.add_axes([0.45,.45,.05,.075])
mr_eq_text = mr_eq_ax.text(0,0,f'$\Omega_m=\Omega_r:\; z={mr_eq_val:.3f}$')
mr_eq_ax.set_axis_off()

m_max_val = getRedshift(NAxis[indexMPeak])
m_max_ax = fig.add_axes([0.45,.425,.05,.075])
m_max_text = m_max_ax.text(0,0,f'max$(\Omega_m):\; z={m_max_val:.3f}$')
m_max_ax.set_axis_off()

msf_eq_val = getRedshift(NAxis[indexMPhi_eq])
msf_eq_ax = fig.add_axes([0.45,.4,.05,.075])
msf_eq_text = msf_eq_ax.text(0,0,f'$\Omega_m=\Omega_\phi:\; z={msf_eq_val:.3f}$')
msf_eq_ax.set_axis_off()



todayLine, = dens_ax.plot([0,0], [-0.2,1.2], 'k', label='Today:\; $z = 0$')
MR_eqLine, = dens_ax.plot([NAxis[indexMR_eq],NAxis[indexMR_eq]], [-0.2,1.2], 'k:', linewidth = 0.75,
                          label = f'$\Omega_m=\Omega_r:\; z={mr_eq_val:.1f}$')
MPeakLine, = dens_ax.plot([NAxis[indexMPeak],NAxis[indexMPeak]], [-0.2,1.2], 'k--', linewidth = 0.75,
                          label = f'max$(\Omega_m):\; z={m_max_val:.1f}$')
MPhi_eqLine, = dens_ax.plot([NAxis[indexMPhi_eq],NAxis[indexMPhi_eq]], [-0.2,1.2], 'k-.', linewidth = 0.75,
                          label = f'$\Omega_m=\Omega_\phi:\; z={msf_eq_val:.1f}$')

Radn_dens_plot, = dens_ax.plot(NAxis, rad_dens, 'r')
Mass_dens_plot, = dens_ax.plot(NAxis, mass_dens, 'g')
Phi_dens_plot, = dens_ax.plot(NAxis, phi_dens, 'b')



rScalingLine,  = gamma_ax.plot([N[-1], N[0]], [4/3, 4/3], "r:", linewidth = .75, label='$\gamma_r=4/3$')
mScalingLine,  = gamma_ax.plot([N[-1], N[0]], [1, 1], "g--", linewidth = .75, label='$\gamma_m=1$')
accel_plot, = gamma_ax.plot(NAxis,
                accelerationExpression(pathx,pathy,pathz),'darkorange', label='Acceleration')
effective_eos, = gamma_ax.plot(NAxis, gamma_phi(pathx, pathy), 'b-', linewidth=1, label = r'$\gamma_\phi$')
CCScalingLine, = gamma_ax.plot([N[-1], N[0]], [0, 0], "k-.", linewidth = .75, label='$\gamma_\Lambda=0$')



x_i, y_i, z_i = state_0[0], state_0[1], state_0[2]
state0_point, = track_ax.plot(x_i,y_i,z_i, 'cX')

track = track_ax.plot(
                pathx, pathy, pathz, 'b', linewidth=2)[0]

fixedPoints, fixedPoints_labels = fixedPoints_func(lam)
for point in fixedPoints:
    plot, = track_ax.plot(point[0], point[1], point[2], 'or', markeredgewidth=0.5, markeredgecolor='black')
    fixedPoint_plots.append(plot)

track.set_visible(True)


xRunningIntegral = np.append(0, cumulative_trapezoid(pathx, NAxis))

hubbleFromY = HfromY(pathy, xRunningIntegral, lam)
hubbleFromY *= H_0 / hubbleFromY[indexToday]

rho_r = rad_dens * hubbleFromY**2
rho_m = mass_dens * hubbleFromY**2
rho_phi = phi_dens * hubbleFromY**2

z_decouple = 1100
N_decouple = -np.log(1+z_decouple)


N_start_idx = np.argmax(N > N_plot_min)
N_end_idx   = np.argmax(N > N_plot_max)


hubble_ax.plot([0,0],[1e0,1e20],'k-', linewidth=1, label='Today')
hubble_ax.plot([N_decouple,N_decouple],[1e0,1e20],'darkorange',
                linestyle=':', linewidth=1, label=r'$z_{\mathrm{dec}}\approx 1100$')
hubble_plot, = hubble_ax.plot(NAxis[N_start_idx:N_end_idx],
                               hubbleFromY[N_start_idx:N_end_idx], color='b',label=r'$H^{(\phi)}(N)$')
#hubble_plot, = hubble_ax.plot(NAxis, hubbleFromY1)

hubble_SN_line = hubble_ax.plot([N_plot_min,N_plot_max], [H0_SN, H0_SN],
                                 alpha=0.75, color = "cyan", label=r'$H_0^{\mathrm{SN}}=73.04\pm 1.04$')
hubble_PL_line = hubble_ax.plot([N_plot_min,N_plot_max], [H0_PL, H0_PL],
                                 alpha=0.75, color = "magenta", label=r'$H_0^{\mathrm{PL}}=67.85\pm 0.52$')


#_____________________Setting plot labels etc________________________________

'''GUI Settings'''


track_ax.set(xlabel='$x$', ylabel='$y$', zlabel='$z$',
             xlim = [-1,1], ylim = [0,1], zlim = [0,1],
             xticks = [-1, -0.5, 0, 0.5, 1],
             yticks = [0, 0.5, 1],
             zticks = [0, 0.5, 1])
track_ax.set_box_aspect([2, 1, 1])
track_ax.axis("off")


gamma_ax.set(yticks = [-1, 0, 1, 4/3, 2], ylim=[-1.1,2.25],
            yticklabels = ['$-1$', '$0$','$1$', '$4/3$', '$2$'],
            xlim = [N_plot_min, N_plot_max], xticks=[-15,-10,-5,0,5])
gamma_ax.yaxis.set_ticks_position('both')
gamma_ax.legend(fontsize=12, loc='upper right', ncol=2)
gamma_ax.tick_params(axis='x', which='both', labelbottom=True) 
gamma_ax.set_xlabel('$N$')
gamma_ax.set_ylabel('$\gamma$', rotation = 0)

gam2_ax = gamma_ax.twinx()
gam2_ax.set_ylabel('Acceleration')
gam2_ax.set(yticks = [-1, 0, 1, 4/3, 2], ylim=[-1.1,2.25],
            yticklabels = ['','','','',''])


dens_ax.set(xlabel="$N$", ylabel="Density Parameters",
            ylim=[-0.1,1.2],yticks=[0,1/4,1/2,3/4,1],
            yticklabels = ['$0$','$1/4$','$1/2$', '$3/4$', '$1$'],
            xlim = [N_plot_min, N_plot_max], xticks=[-15,-10,-5,0,5])

hubble_ax.set_yscale('log', base=10, subs=[10**x
                         for x in (0.25, 0.5, 0.75)], nonpositive='mask')
hubble_ax.set(  xlabel = "$N$",
                ylabel=r"$\log_{10}(H)$",
                xlim = [N_plot_min, N_plot_max],
                ylim = [1e0,1e15],
                yticks=[1e0,1e5,1e10,1e15],
                yticklabels=['$0$', '$5$', '$10$','$15$'],
                xticks=[-15,-10,-5,0,5])
hubble_ax.legend(loc='upper right')


#Run the code
window_gui.mainloop()
#End
