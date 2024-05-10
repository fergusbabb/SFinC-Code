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
plt.rcParams.update({"text.usetex": True, "font.family": "serif",
                     "font.serif": ["Computer Modern Serif"]})
plt.rc('axes', labelsize=14, titlesize=15)
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12



#_________________________Set up main window________________________________
#Main GUI figure
window = tk.Tk()
window.title('GUI for Matter and Radiation')
window.geometry('1600x950')
fig = Figure(figsize=(16, 9.5)) #1600x950 pixels
fig.set_facecolor('white')

#Define Axes
#Tracks plot figure
track_axis_dims = [-.025,.525,.45,.5]
track_ax = fig.add_axes(track_axis_dims, projection='3d')
track_ax.view_init(elev=24, azim=66)

#Colour bar axes
cbar_ax_dims = [.375,.6,.015,.35]
cbar_ax = fig.add_axes(cbar_ax_dims)

#Relative Density axes
dens_axis_dims = [.625,.125,.35,.275]
dens_ax = fig.add_axes(dens_axis_dims)

#Acceleration axes
accel_axis_dims = [.625,.4,.35,.275]
accel_ax = fig.add_axes(accel_axis_dims)

#EoS Axes
gamma_axis_dims = [.625,.675,.35,.275]
gamma_ax = fig.add_axes(gamma_axis_dims)

#Hubble plot Axes
d_lum_ax_dims = [.0675,.125,.35,.3675] 
d_lum_ax = fig.add_axes(d_lum_ax_dims)


'''
To produce report figures comment above axes, uncomment below
'''

''' You cant loop this for some reason....'''

# window_4_report = tk.Tk()
# window_4_report.title('Window to generate plots for report 1')
# window_4_report.geometry('750x500')
# fig2 = Figure(figsize=(7.5, 5)) #750x500 pixels
# fig2.set_facecolor('white')

# window_4_report2 = tk.Tk()
# window_4_report2.title('Window to generate plots for report 2')
# window_4_report2.geometry('750x500')
# fig3 = Figure(figsize=(7.5, 5)) #750x500 pixels
# fig3.set_facecolor('white')

# window_4_report3 = tk.Tk()
# window_4_report3.title('Window to generate plots for report 3')
# window_4_report3.geometry('750x500')
# fig4 = Figure(figsize=(7.5, 5)) #750x500 pixels
# fig4.set_facecolor('white')

# window_4_report4 = tk.Tk()
# window_4_report4.title('Window to generate plots for report 4')
# window_4_report4.geometry('750x500')
# fig5 = Figure(figsize=(7.5, 5)) #750x500 pixels
# fig5.set_facecolor('white')

# window_4_report5 = tk.Tk()
# window_4_report5.title('Window to generate plots for report 5')
# window_4_report5.geometry('750x500')
# fig6 = Figure(figsize=(7.5, 5)) #750x500 pixels
# fig6.set_facecolor('white')


# track_axis_dims2 = [0,.075,.9,.9]
# track_ax = fig2.add_axes(track_axis_dims2, projection='3d')
# track_ax.view_init(elev=24, azim=66)

# cbar_ax_dims2 = [.7,.25,.02,.6]
# cbar_ax = fig2.add_axes(cbar_ax_dims2)

# dens_axis_dims2 = [.15,.25,.8,.7]
# dens_ax = fig3.add_axes(dens_axis_dims2)

# accel_axis_dims2 = [.15,.25,.8,.7]
# accel_ax = fig4.add_axes(accel_axis_dims2)

# gamma_axis_dims2 = [.15,.25,.8,.7]
# gamma_ax = fig5.add_axes(gamma_axis_dims2)

# d_lum_ax_dims2 = [.15,.25,.8,.7]
# d_lum_ax = fig6.add_axes(d_lum_ax_dims2)

# windows = [window_4_report, window_4_report2, window_4_report3, window_4_report4, window_4_report5]
# figures = [fig2, fig3, fig4, fig5, fig6]

# def setup_canvas_and_toolbar(figval, parent_window, toolbar_parent=None):
#     if toolbar_parent is None:
#         toolbar_parent = parent_window
#     canvas = FigureCanvasTkAgg(figval, parent_window)
#     canvas.draw()
#     canvas.get_tk_widget().place(relheight=1, relwidth=1)
#     NavigationToolbar2Tk(canvas, toolbar_parent)

# # Apply the setup to each figure and window pair
# for figval, win in zip(figures, windows):
#     setup_canvas_and_toolbar(figval, win)


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

#__________________________Initial values_____________________________
pathnum = 1

lam1_0 = 2
lam1_min = -10
lam1_max = 10

lam2_0 = -10
lam2_min = -10
lam2_max = 10

Ni = -8
NiForward = 16

N = np.linspace(0, Ni, 8000)
NForward = np.linspace(0, NiForward, abs(int(NiForward*1000)))

z = np.exp(-N) - 1
c = 1 #3e5     # Given in km/s
V = z * c   # ""
h = 0.738
H_0 = 100*h # km/(s Mpc)
xinit = np.linspace(-0.99, 0.99, pathnum)

#Observed Values
Omega_phi_0 = 0.738
gamma_phi_value = 0.013
Omega_r0 = 4.984e-5

peakOmega_rConsistent = np.array([])

x0Squared = gamma_phi_value / (2 * Omega_phi_0)
y0Squared = Omega_phi_0 - x0Squared

Omega_m0 = 1 - x0Squared - y0Squared - Omega_r0

#state_0 = [round(np.sqrt(x0Squared),3),
#           round(np.sqrt(y0Squared),3),
#           round(np.sqrt(y0Squared) * 1e-1,3),
#           round(np.sqrt(Omega_r0),3)]

state_0 = [0.0000001,
           0.00001,
           0.00001,
           0.975]


#_____________________________Define ODE for 2 fluids_______________________

def ODEs(state, N, lam):
    # coords[0]    corresponds to the variable x
    # coords[1:-1] corresponds to the variables yn
    # coords[-1]   corresponds to the variable z

    x = state[0]
    y = np.array(state[1:-1])
    z = state[-1]

    HdotSection = 1 + x**2 - np.sum(y**2) + (1/3)*z**2

    xprime = (-3*x) + (np.sqrt(1.5) *  np.sum(lam * y**2)) + (1.5*x * HdotSection)
    yprime =          (np.sqrt(1.5) *        -lam * x * y) + (1.5*y * HdotSection)
    zprime = (-2*z)                                        + (1.5*z * HdotSection)

    return np.append(xprime, np.append(yprime, zprime))

#___________________________Function for fixed points________________________
def fixedPoints_func(lam):
    fixedPoints_labels = ['$O$', '$A^+$', '$A^-$']

    fixedPoints = np.array([
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0]
    ])
    
    for n, lambdaValue in enumerate(lam):
        # B1 and B2
        if lambdaValue**2 < 6:
            fixedPoints = np.append(fixedPoints,
                            [[lambdaValue/np.sqrt(6),
                                np.sqrt(1 - (lambdaValue**2)/6), 0]], axis=0)

            fixedPoints_labels.append('$B' + str(n+1) + '$')

        # C1 and C2
        if lambdaValue**2 > 3: #Makes sure to only plot positive y points
            fixedPoints = np.append(fixedPoints,
                                    [[np.sqrt(3/2)/lambdaValue,
                                        np.sqrt(3/2)/lambdaValue, 0]], axis=0)
            fixedPoints_labels.append('$C' + str(n+1) + '$')

        # E1 and E2
        if lambdaValue >= 0: #Makes sure to only plot positive y points
            if lambdaValue**2 > 4:
                fixedPoints = np.append(fixedPoints,
                                    [[2 * np.sqrt(2/3) / lambdaValue,
                                    2 / (lambdaValue * np.sqrt(3)),
                                    np.sqrt(1 - (4/lambdaValue**2))]], axis=0)
                fixedPoints_labels.append('$E' + str(n+1) + '$')

    fixedPoints = np.append(fixedPoints, [[0, 0, 1]], axis=0)  
    fixedPoints_labels.append('$D$')
    fixedPoints = np.append(fixedPoints, [[0, 1, 0]], axis=0)  
    fixedPoints_labels.append('$F$')  
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


#_______________________________Update track plots___________________________
def update_plot(event):
    #Before update_plot is called lam sliders are updated
    #Sliders will always have current lambda/ gamma values
    lam1 = lambda1_slide.get()
    lam2 = lambda2_slide.get()

    state_0 = [float(x_entry_val.get()), 
               float(y1_entry_val.get()),
               float(y2_entry_val.get()), 
               float(z_entry_val.get())]
    
    for plot in fixedPoint_plots:
        plot.remove()
    fixedPoint_plots.clear()  # Clear the list of plot objects

    # Get new fixed points and plot them
    fixedPoints, fixedPoints_labels = fixedPoints_func([lam1,lam2])
    for i, point in enumerate(fixedPoints):
        plot, = track_ax.plot(point[0], point[1], point[2], 'or')
        fixedPoint_plots.append(plot)
    
    x_i, y_i, z_i = state_0[0], np.sqrt(state_0[1]**2 + state_0[2]**2), state_0[3]
    #state0_point.set_data(x_i, y_i)
    #state0_point.set_3d_properties(z_i)
    
    #Plot all paths with updated values


    #Update the quiver vectors
    quiver_vectors = np.array([ODEs([pt[0], pt[1], pt[2]], N, lam1)
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
    solution = odeint(ODEs, state_0, N, args = (np.array([lam1, lam2]),))
    pathx  = solution[:, 0]
    pathy1 = solution[:, 1]
    pathy2 = solution[:, 2]
    pathz  = solution[:, 3]
    
    #Add both ways
    solutionForward = odeint(ODEs, state_0, NForward, args = (np.array([lam1, lam2]),))
    pathx  = np.append(pathx[::-1],  solutionForward[:, 0])
    pathy1 = np.append(pathy1[::-1], solutionForward[:, 1])
    pathy2 = np.append(pathy2[::-1], solutionForward[:, 2])
    pathz  = np.append(pathz[::-1],  solutionForward[:, 3])
    pathy = np.sqrt(pathy1**2 + pathy2**2)

    #Update tracks
    track.set_data(pathx, pathy)
    track.set_3d_properties(pathz)

    #Update acceleration plot
    accel_plot_new_data = accelerationExpression(
                                pathx,pathy,pathz)
    accel_plot.set_ydata(accel_plot_new_data)

    
    #Find when in N each event occured
    NAxis = np.append(N[::-1], NForward)
    indexToday = np.argmin(np.abs(phi_dens-Omega_phi_0))
    

    indexMR_eq = np.argmin(np.abs(mass_dens-rad_dens)[0:indexToday])
    indexMPhi_eq = np.argmin(np.abs(mass_dens-phi_dens)[indexMR_eq:-1]) + indexMR_eq
    indexMPeak = np.argmax(mass_dens)
    rScalingLine.set_xdata([NAxis[0], NAxis[-1]])
    mScalingLine.set_xdata([NAxis[0], NAxis[-1]])
    #gamma_ax.set(xlim=[NAxis[0], NAxis[-1]])
    NAxis -= NAxis[indexToday]

    MR_eqLine.set_xdata([NAxis[indexMR_eq], NAxis[indexMR_eq]])
    MPhi_eqLine.set_xdata([NAxis[indexMPhi_eq], NAxis[indexMPhi_eq]])
    MPeakLine.set_xdata([NAxis[indexMPeak], NAxis[indexMPeak]])

    #Update relative density plots
    Radn_dens_plot.set_ydata(pathz**2)
    Mass_dens_plot.set_ydata(1 - pathx**2 - pathy**2 - pathz**2)
    Phi_dens_plot.set_ydata(pathx**2 + pathy**2)
    #y1_dens_plot.set_ydata(pathy1**2)
    #y2_dens_plot.set_ydata(pathy2**2)


    mr_eq_val = getRedshift(NAxis[indexMR_eq])
    mr_eq_text.set_text(f'$\Omega_m=\Omega_r:\; z={mr_eq_val:.3f}$')

    m_max_val = getRedshift(NAxis[indexMPeak])
    m_max_text.set_text(f'max$(\Omega_m):\; z={m_max_val:.3f}$')

    msf_eq_val = getRedshift(NAxis[indexMPhi_eq])
    msf_eq_text.set_text(f'$\Omega_m=\Omega_\phi:\; z={msf_eq_val:.3f}$')
    #rsf_eq_val = 
    #rsf_eq_text.set_text(f'$\Omega_\phi=\Omega_r:\; z={rsf_eq_val:.3f}$')


    #Update EoS plot
    path_gamma_phi = gamma_phi(pathx, pathy) 
    effective_eos.set_ydata(path_gamma_phi)
    
    #Update redshift plot
    d_L = (c) * (1 + z) * odeint(
    d_L_IntegrandScalar, 0, z, args=(
        zAxis, Omega_m0, Omega_r0, Omega_phi_0, path_gamma_phi
        )).transpose()[0]
    
    integral_plot.set_ydata(d_L)

    #Show plots
    fig.canvas.draw()



#____________________________Defining Widgets________________________________
#When cursor clicks on region plot update to clicked value

#Canvas is where figure is placed to window
canvas = FigureCanvasTkAgg(fig, window)
canvas.draw() #Show canvas (ie show figure)

#Show the navigation toolbar
NavigationToolbar2Tk(canvas, window)


#Lambda Slider Labels. Initialise, place, and hide weird borders
lambda1_label_ax = fig.add_axes([.455,.9525,.05,.05])
lambda1_label_ax.text(0,0,'$\lambda_1$ value')
lambda1_label_ax.set_axis_off()

lambda2_label_ax = fig.add_axes([.515,.9525,.05,.05])
lambda2_label_ax.text(0,0,'$\lambda_2$ value')
lambda2_label_ax.set_axis_off()

#Lambda Slider. Initialise, add interaction, place, hide borders
lambda1_slide = tk.Scale(window, from_ = lam1_min, to = lam1_max, width = 20, length = 250, resolution=0.001)
lambda1_slide.set(lam1_0)
lambda1_slide.bind("<ButtonRelease-1>", update_plot)
lambda1_slide.place(relx=0.45, rely=0.05, relheight=0.35, relwidth=0.05)
lambda1_slide.configure(bg = 'white', borderwidth=0)

lambda2_slide = tk.Scale(window, from_ = lam1_min, to = lam1_max, width = 20, length = 250, resolution=0.001)
lambda2_slide.set(lam2_0)
lambda2_slide.bind("<ButtonRelease-1>", update_plot)
lambda2_slide.place(relx=0.51, rely=0.05, relheight=0.35, relwidth=0.05)
lambda2_slide.configure(bg = 'white', borderwidth=0)

def submit():
    x0 = float(x_entry_val.get())
    y01 = float(y1_entry_val.get())
    y02 = float(y2_entry_val.get())
    z0 = float(z_entry_val.get())

    if x0**2 + y01**2 + y02**2 + z0**2 <= 1:
        update_plot(0)
    else:
        print('Not valid values, check $x^2 + y1^2 + y02^2 + z^2 <= 1$')

x_entry_val=tk.StringVar()
y1_entry_val=tk.StringVar()
y2_entry_val=tk.StringVar()
z_entry_val=tk.StringVar()

x_entry_label_ax = fig.add_axes([0.045,.525,.05,.075])
x_entry_label_ax.text(0,0,'$x_0$:')
x_entry_label_ax.set_axis_off()

y1_entry_label_ax = fig.add_axes([0.1225,.525,.05,.075])
y1_entry_label_ax.text(0,0,'$y_{01}$:')
y1_entry_label_ax.set_axis_off()

y2_entry_label_ax = fig.add_axes([0.2025,.525,.05,.075])
y2_entry_label_ax.text(0,0,'$y_{02}$:')
y2_entry_label_ax.set_axis_off()

z_entry_label_ax = fig.add_axes([0.285,.525,.05,.075])
z_entry_label_ax.text(0,0,'$z_0$:')
z_entry_label_ax.set_axis_off()


x_entry = tk.Entry(window, textvariable = x_entry_val, relief='solid'
            ).place(relx=0.06, rely=0.457, relheight=0.03, relwidth=0.05)
y1_entry = tk.Entry(window, textvariable = y1_entry_val, relief='solid'
            ).place(relx=0.14, rely=0.457, relheight=0.03, relwidth=0.05)
y2_entry = tk.Entry(window, textvariable = y2_entry_val, relief='solid'
            ).place(relx=0.22, rely=0.457, relheight=0.03, relwidth=0.05)
z_entry = tk.Entry(window, textvariable = z_entry_val, relief='solid'
            ).place(relx=0.3, rely=0.457, relheight=0.03, relwidth=0.05)

x_entry_val.set(state_0[0])    
y1_entry_val.set(state_0[1])
y2_entry_val.set(state_0[2])
z_entry_val.set(state_0[3])

sub_btn=tk.Button(window, text = 'Submit', command = submit
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
quiver_vectors = np.array([ODEs([pt[0], pt[1], pt[2]], N, lam1_0)
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

#_______________________________Hubble Fill Regions___________________________
w_Lam_0 = -1
w_pos_err = 0.15
w_neg_err = -0.15

def plot_d_luminosity(ax, z, d_L, d_L_bounds, label, color, fill_alpha=0.2):
    ax.plot(z, d_L, label=label, color=color, lw=2)
    if d_L_bounds is not None:
        ax.fill_between(z, d_L_bounds[0], d_L_bounds[1], alpha=fill_alpha, color=color)
        ax.plot(z, d_L_bounds[0], color=color, lw=0.5, alpha = 0.6)
        ax.plot(z, d_L_bounds[1], color=color, lw=0.5, alpha = 0.6)

def setup_luminosity_plots():
    d_L_for_fill = []
    #Define colors and Omega_Lambda0 values
    configurations = [
        (0.68, 'cyan'),
        (0.73, 'magenta')
    ]

    #Gather all d_L values for bounds and normal plotting
    for Omega_Lambda0, color in configurations:
        temp_list = []
        for w_Lambda in [w_Lam_0 + w_neg_err,
                         w_Lam_0,
                         w_Lam_0 + w_pos_err]:
            d_L = (c) * (1 + z) * odeint(
                d_L_Dark_Energy, 0, z, args=(
                    zAxis, 1-Omega_Lambda0-Omega_r0, Omega_r0, Omega_Lambda0, w_Lambda
                )
            ).transpose()[0]
            temp_list.append(d_L)
        d_L_for_fill.append(temp_list)  #Store d_L values for each configuration
    
    #Plot and fill between using stored d_L values
    for (Omega_Lambda0, color,), d_L_values in zip(configurations, d_L_for_fill):
        # Plot the middle value normally and fill between the bounds
        plot_d_luminosity(d_lum_ax, z, d_L_values[1], [d_L_values[0], d_L_values[2]], 
                          f"$\Omega_{{\Lambda}}^{{(0)}}={Omega_Lambda0},\; w_{{\Lambda}}={w_Lam_0}$", color)


#___________________________________Initial Plot_____________________________

#Initial plot
fixedPoint_plots = []

rScalingLine, = gamma_ax.plot([N[-1], NForward[-1]], [4/3, 4/3], "k-.", linewidth = 0.5, label='$\gamma_r$')

mScalingLine, = gamma_ax.plot([N[-1], NForward[-1]], [1, 1], "k--", linewidth = 0.5, label='$\gamma_m$')

lam1 = lambda1_slide.get()
lam2 = lambda2_slide.get()

# Solve the system of ODEs using odeint
solution = odeint(ODEs, state_0, N, args = (np.array([lam1, lam2]),))
pathx  = solution[:, 0]
pathy1 = solution[:, 1]
pathy2 = solution[:, 2]
pathz  = solution[:, 3]

solutionForward = odeint(ODEs, state_0, NForward, args = (np.array([lam1, lam2]),))
pathx  = np.append(pathx[::-1],  solutionForward[:, 0])
pathy1 = np.append(pathy1[::-1], solutionForward[:, 1])
pathy2 = np.append(pathy2[::-1], solutionForward[:, 2])
pathz  = np.append(pathz[::-1],  solutionForward[:, 3])
pathySquared = pathy1**2 + pathy2**2
pathy = np.sqrt(pathySquared)

path_gamma_phi = gamma_phi(pathx, pathy)

NAxis = np.append(N[::-1], NForward)    

rad_dens = pathz**2
mass_dens = 1 - pathx**2 - pathySquared - pathz**2
phi_dens = pathx**2 + pathySquared
total_dens = rad_dens + mass_dens + phi_dens

indexToday = np.argmin(np.abs(phi_dens-Omega_phi_0))

#Find when in N each event occured
indexMR_eq = np.argmin(np.abs(mass_dens-rad_dens)[0:indexToday])
indexMPhi_eq = np.argmin(np.abs(mass_dens-phi_dens)[indexMR_eq:-1]) + indexMR_eq
indexMPeak = np.argmax(mass_dens)
NAxis -= NAxis[indexToday]
zAxis = getRedshift(NAxis)


#rScalingLine.set_xdata([NAxis[0], NAxis[-1]])
#mScalingLine.set_xdata([NAxis[0], NAxis[-1]])
#gamma_ax.set(xlim=[NAxis[0], NAxis[-1]])



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

#rsf_eq_val = 
#rsf_eq_ax = fig.add_axes([0.45,.375,.05,.075])
#rsf_eq_text = rsf_eq_ax.text(0,0,f'$\Omega_\phi=\Omega_r:\; z={rsf_eq_val:.3f}$')
#rsf_eq_ax.set_axis_off()

todayLine, = dens_ax.plot([0,0], [-0.2,1.2], 'k', label='Today:\; $z = 0$')
MR_eqLine, = dens_ax.plot([NAxis[indexMR_eq],NAxis[indexMR_eq]], [-0.2,1.2], 'k:', linewidth = 0.75,
                          label = f'$\Omega_m=\Omega_r:\; z={mr_eq_val:.3f}$')
MPeakLine, = dens_ax.plot([NAxis[indexMPeak],NAxis[indexMPeak]], [-0.2,1.2], 'k--', linewidth = 0.75,
                          label = f'max$(\Omega_m):\; z={m_max_val:.3f}$')
MPhi_eqLine, = dens_ax.plot([NAxis[indexMPhi_eq],NAxis[indexMPhi_eq]], [-0.2,1.2], 'k-.', linewidth = 0.75,
                          label = f'$\Omega_m=\Omega_\phi:\; z={msf_eq_val:.3f}$')

Radn_dens_plot, = dens_ax.plot(NAxis, rad_dens, 'r',
        label = "$\Omega_r$")
Mass_dens_plot, = dens_ax.plot(NAxis, mass_dens, 'g',
        label = "$\Omega_m$")
Phi_dens_plot, = dens_ax.plot(NAxis, phi_dens, 'b',
        label = "$\Omega_\phi$")

#y1_dens_plot, =  dens_ax.plot(NAxis, pathy1**2, 'b--',
#        label = "$y_1^2$")
#y2_dens_plot, =  dens_ax.plot(NAxis, pathy2**2, 'b--',
#        label = "$y_2^2$")



x_i, y_i, z_i = state_0[0], np.sqrt(state_0[1]**2 + state_0[2]**2), state_0[3]


# state0_point, = track_ax.plot(x_i,y_i,z_i, 'cX')
track = track_ax.plot(
                pathx, pathy, pathz, 'b', linewidth=2)[0]
accel_plot, = accel_ax.plot(NAxis,
                accelerationExpression(pathx,pathy,pathz),'orange')
effective_eos, = gamma_ax.plot(NAxis, gamma_phi(pathx, pathy), 'springgreen',
        label = r'$\gamma_\phi = 2x^2/(x^2+y^2)$')

fixedPoints, fixedPoints_labels = fixedPoints_func([lam1,lam2])
for point in fixedPoints:
    plot, = track_ax.plot(point[0], point[1], point[2], 'or')
    fixedPoint_plots.append(plot)


track.set_visible(True)


d_L = (c) * (1 + z) * odeint(
    d_L_IntegrandScalar, 0, z, args=(
        zAxis, Omega_m0, Omega_r0, Omega_phi_0, path_gamma_phi
        )).transpose()[0]


#Call hubble fill function before changing plot, so it is below
setup_luminosity_plots()

integral_plot, = d_lum_ax.plot(z, d_L,
                    label = f"$\Omega_{{\phi}}^{{(0)}} = {Omega_phi_0}$", color = 'b', linewidth=2)



#_____________________Setting plot labels etc________________________________

'''GUI Settings'''


track_ax.set(xlabel='$x$', ylabel='$y$', zlabel='$z$',
             xlim = [-1,1], ylim = [0,1], zlim = [0,1],
             xticks = [-1, -0.5, 0, 0.5, 1],
             yticks = [0, 0.5, 1],
             zticks = [0, 0.5, 1])
track_ax.set_box_aspect([2, 1, 1])
track_ax.axis("off")


accel_ax.set(ylabel="Acceleration", ylim=[-1.1,1.1],
             yticks=[-1,-1/2,0,1/2,1], yticklabels = ['$-1$','$-1/2$', '$0$', '$1/2$', '$1$'],
             xlim=[-8,3], xticks = [-8,-6,-4,-2,0,2],
             xticklabels = ['$-8$', '$-6$', '$-4$', '$-2$','$0$','$2$'])
accel_ax.tick_params(axis='x', which='both', labelbottom=False)
accel_ax.yaxis.set_ticks_position('both')
static_line = accel_ax.plot([-8,3],[0,0], "k--", linewidth = 0.5)


gamma_ax.set(ylabel="$\gamma_\phi$", yticks = [0, 1, 4/3, 2], ylim=[-0.1,2.1],
             yticklabels = ['$0$','$1$', '$4/3$', '$2$'], xlim=[-8,3], xticks = [-8,-6,-4,-2,0,2],
             xticklabels = ['$-8$', '$-6$', '$-4$', '$-2$','$0$','$2$'])
gamma_ax.tick_params(axis='x', which='both', labelbottom=False)
gamma_ax.yaxis.set_ticks_position('both')

dens_ax.set(xlabel="$N$", ylabel="Density Parameters",
             ylim=[-0.1,1.1],yticks=[0,1/4,1/2,3/4,1],
             yticklabels = ['$0$','$1/4$','$1/2$', '$3/4$', '$1$'],
             xlim=[-8,3], xticks = [-8,-6,-4,-2,0,2],
             xticklabels = ['$-8$', '$-6$', '$-4$', '$-2$','$0$','$2$'])

d_lum_ax.set(ylabel = "$H_0d_L$ ", xlabel= '$z$',
              xlim=[0,3], ylim=[0,6],
              xticks=[0,1,2,3], yticks=[0,1,2,3,4,5,6],
              xticklabels = ['$0$','$1$','$2$', '$3$'],
              yticklabels = ['$0$','$1$','$2$', '$3$', '$4$','$5$','$6$'])



'''Figure Production Settings'''
# track_ax.set(xlabel='$x$', ylabel='$y$', zlabel='$z$',
#              xlim = [-1,1], ylim = [0,1], zlim = [0,1],
#              xticks = [-1, -0.5, 0, 0.5, 1],
#              yticks = [0, 0.5, 1],
#              zticks = [0, 0.5, 1])
# track_ax.set_box_aspect([2, 1, 1])
# track_ax.axis("off")



# static_line = accel_ax.plot([-8,3],[0,0], "k--", linewidth = 0.5)
# accel_ax.set(ylabel="Acceleration", ylim=[-1.1,1.1],
#              yticks=[-1,-1/2,0,1/2,1], yticklabels = ['$-1$','$-1/2$', '$0$', '$1/2$', '$1$'],
#                xlim=[-8,3], xticks = [-8,-6,-4,-2,0,2], xlabel="$N$",
#             xticklabels = ['$-8$', '$-6$', '$-4$', '$-2$','$0$','$2$'])
# accel_ax.yaxis.set_ticks_position('both')
# accel_ax.tick_params(axis='x', which='both', labelbottom=True)


# gamma_ax.set(ylabel="$\gamma_\phi$", yticks = [0, 1, 4/3, 2], ylim=[-0.1,2.1], xlabel="$N$",
#             yticklabels = ['$0$','$1$', '$4/3$', '$2$'], xlim=[-8,3], xticks = [-8,-6,-4,-2,0,2],
#             xticklabels = ['$-8$', '$-6$', '$-4$', '$-2$','$0$','$2$'])
# gamma_ax.yaxis.set_ticks_position('both')
# gamma_ax.legend(fontsize=12, loc='best')
# gamma_ax.tick_params(axis='x', which='both', labelbottom=True) 


# dens_ax.set(xlabel="$N$", ylabel="Density Parameters",
#             ylim=[-0.1,1.1],yticks=[0,1/4,1/2,3/4,1],
#             yticklabels = ['$0$','$1/4$','$1/2$', '$3/4$', '$1$'],
#             xlim=[-8,3], xticks = [-8,-6,-4,-2,0,2],
#             xticklabels = ['$-8$', '$-6$', '$-4$', '$-2$','$0$','$2$'])

# #Additional code for making paper plots
# legend_lines1 = []
# legend_lines1.append([todayLine, MR_eqLine, MPeakLine, MPhi_eqLine])
# legend_lines2 = []
# legend_lines2.append([Radn_dens_plot, Mass_dens_plot, Phi_dens_plot])

# dens_ax.text(NAxis[indexMR_eq]-0.5,.1,f'$z={mr_eq_val:.1f}$',backgroundcolor='1', fontsize=12)
# dens_ax.text(NAxis[indexMPeak]-1,1,f'$z={m_max_val:.1f}$',backgroundcolor='1', fontsize=12)
# dens_ax.text(NAxis[indexMPhi_eq]-0.5,1,f'$z={msf_eq_val:.1f}$',backgroundcolor='1', fontsize=12)

# legend1 = dens_ax.legend(legend_lines1[0], ["Today","$\Omega_m=\Omega_r$","max$(\Omega_m)$",
#                                             "$\Omega_m=\Omega_\phi$"], loc='upper left', fontsize=12)
# legend2 = dens_ax.legend(legend_lines2[0], ['$\Omega_r$', '$\Omega_m$', '$\Omega_\phi$'],
#                          loc='center left', bbox_to_anchor=(0, .45), fontsize=12)
# dens_ax.add_artist(legend1)
# dens_ax.yaxis.set_ticks_position('both')



# d_lum_ax.set(ylabel = "$H_0d_L$", xlabel= '$z$',
#               xlim=[0,3], ylim=[0,6],
#               xticks=[0,1,2,3], yticks=[0,1,2,3,4,5,6],
#               xticklabels = ['$0$','$1$','$2$', '$3$'],
#               yticklabels = ['$0$','$1$','$2$', '$3$', '$4$','$5$','$6$'])
# d_lum_ax.legend(loc=4)



# fig2.savefig("Figures/One Potential/ Track_lambda2_{}.svg".format(round(lam_0**2)), format='svg')
# fig3.savefig("Figures/One Potential/ Density_lambda2_{}.svg".format(round(lam_0**2)), format='svg')
# fig4.savefig("Figures/One Potential/ Accel_lambda2_{}.svg".format(round(lam_0**2)), format='svg')
# fig5.savefig("Figures/One Potential/ Gamma_lambda2_{}.svg".format(round(lam_0**2)), format='svg')
# fig6.savefig("Figures/One Potential/ Hubble_lambda2_{}.svg".format(round(lam_0**2)), format='svg')



#Run the code
window.mainloop()
