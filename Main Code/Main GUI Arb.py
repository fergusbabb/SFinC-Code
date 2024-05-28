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
plt.rc('axes', labelsize=14, titlesize=15)
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

#_________________________Set up main window___________________________________
#Standard tkinter window set up
window = tk.Tk()
window.title('GUI for Arbitrary Fluid')
window.geometry('1200x800')



#Tracks plot figure
fig = Figure(figsize=(12, 8)) #1200x900 pixels
fig.set_facecolor('white')
track_axis_dims = [.075,.6,.35,.3]
track_ax = fig.add_axes(track_axis_dims, aspect='equal')

#Colour bar axes
cbar_ax_dims = [.45,.6,.015,.3]
cbar_ax = fig.add_axes(cbar_ax_dims)

gamma_ax_dims = [.6,.6,.35,.3]
gamma_ax = fig.add_axes(gamma_ax_dims)

#Hubble plot Axes
lam_gam_dims = [.225,.15,.45,.35] 
lam_gam_ax = fig.add_axes(lam_gam_dims)



#Bounding Circle
theta = np.linspace(0, np.pi, 150)
track_ax.plot(np.cos(theta), np.sin(theta),
         'k', linewidth=1)
track_ax.plot([-1,1], [0,0], 'k', linewidth=1)

canvas = FigureCanvasTkAgg(fig, window) 
canvas.draw() 

#Place Canvas
canvas.get_tk_widget().place(relheight=1,relwidth=1)

#Initial values
gam_0 = 1
lam_0 = 3
pathnum = 50
N = np.linspace(0, 15, 10000) 
xinit = np.linspace(-0.99, 0.99, pathnum)


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

#___________________________Function for fixed points________________________
def fixedPoints_func(lam, gam):
    fixedPoints_labels = ['$O$', '$A^+$', '$A^-$','$B$']

    fixedPoints = np.array([
    [0, 0],
    [1, 0],
    [-1, 0]
    ])

    if lam**2 < 6:
        fixedPoints = np.append(fixedPoints,
                                [[lam/np.sqrt(6),
                                    np.sqrt(1 - (lam**2)/6)]], axis=0)
        fixedPoints_labels.append('$C$')

    if lam**2 > 3*gam:
        fixedPoints = np.append(fixedPoints,
                                [[np.sqrt(3/2)*gam/lam,
                                    np.sqrt(3*(2-gam)*gam/2)/lam]], axis=0)
        fixedPoints_labels.append('$D$')
    return fixedPoints, fixedPoints_labels

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
    return track_ax.fill_between(x, boundary, np.sqrt(1-x**2),
                                 where = xWhere, alpha=0.2, color = 'orange', label = 'Accelerating region')
    

#_______________________________Effective Eos Parameter______________________
def gamma_phi(x, y):
    return (2*x**2) / (x**2 + y**2)
 
#_______________________________Update track plots___________________________
def update_plot(event):
    #Before update_plot is called lam and gam sliders are updated
    #Sliders will always have current lambda/ gamma values
    lam = lambda_slide.get()
    gam =  gamma_slide.get()

    global fill
    fill.remove()
    fill = compute_fill(gam)# = plt.fill_between(x, acceleration_y_boundary(x), np.sqrt(1-x**2), where = xWhere)

    for plot in fixedPoint_plots:
        plot.remove()
    fixedPoint_plots.clear()  # Clear the list of plot objects

    # Get new fixed points and plot them
    fixedPoints, fixedPoints_labels = fixedPoints_func(lam, gam)
    for i, point in enumerate(fixedPoints):
        plot, = track_ax.plot(point[0], point[1], 'or', markeredgewidth=0.5, markeredgecolor='black')
        fixedPoint_plots.append(plot)
    
    #Update cursor star point
    clickpoint.set_data(lam**2, gam)
    
    #Plot all paths with updated values
    for i in range(0, pathnum):
        initial_conditions = [xinit[i], yinit[i]]  #Initial values for x and y
        
        #Solve the system of ODEs using odeint
        solution = odeint(ODEs, initial_conditions, N, args = (lam, gam))

        solution_x = solution[:, 0]
        solution_y = solution[:, 1]

        main_tracks[i].set_data(solution_x, solution_y)
        if i == 1:
            eos_plots[i].set_ydata(gamma_phi(solution_x, solution_y))
        if i % 1 == 0:
            eos_plots[i].set_ydata(gamma_phi(solution_x, solution_y))

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


    fig.suptitle(f'$\lambda^2 = {lambda_slide.get()**2:.0f},\; \gamma={gamma_slide.get():.0f}$',fontsize=14)

    gam_ScaleLine.set_ydata([gam,gam])
    
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
#This function ONLY done with Chat GPT

# Parameters
radius = 0.95
num_points = 400  # Number of points on the semicircle

# Generate points on a semicircle using Fibonacci-like spacing
points = fibonacci_semicircle(radius, num_points)
points = np.vstack((points,[0,0.01]))


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
clickpoint, = lam_gam_ax.plot(lam_0**2, gam_0, 'r*', alpha=0)

#Define function for interaction with plot
def regions_plot(event):
    if event.inaxes == lam_gam_ax:
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

#Lambda Slider. Initialise, add interaction, place, hide borders

lambda_label_ax = fig.add_axes([.765,.505,.05,.05])
lambda_label_ax.text(0,0,'$\lambda$ value')
lambda_label_ax.set_axis_off()

lambda_slide = tk.Scale(window, from_ = 0.01, to = np.sqrt(12),
                       width = 20, length = 250, resolution=0.001)
lambda_slide.set(lam_0)
lambda_slide.bind("<ButtonRelease-1>", update_plot)
lambda_slide.place(relx=0.75, rely=0.5, relheight=0.35, relwidth=0.075)
lambda_slide.configure(bg = 'white', borderwidth=0)

#Gamma slider. Initialise, add interaction, place, hide borders
gamma_label_ax = fig.add_axes([.865,.505,.05,.05])
gamma_label_ax.text(0,0,'$\gamma$ value')
gamma_label_ax.set_axis_off()

gamma_slide = tk.Scale(window, from_ = 0.01, to = 2,
                       width = 20, length = 250, resolution=0.001)
gamma_slide.set(gam_0)
gamma_slide.bind("<ButtonRelease-1>", update_plot)
gamma_slide.place(relx=0.85, rely=0.5, relheight=0.35, relwidth=0.075)
gamma_slide.configure(bg = 'white', borderwidth=0)

#___________________________________Initial Plot____________________________________

#Define log a linspace for ode
n = np.linspace(0, 4, 5000)

def initial_points(radius, semicirc_num, diam_num):
    # Correct the use of num_points which is not defined
    delta_theta = np.pi / (semicirc_num - 1)
    points = []
    for i in range(semicirc_num):
        theta = i * delta_theta  # Angle for the point
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        points.append([x, y])  # Use list instead of tuple for easier manipulation

    points = np.array(points)  # Convert to numpy array with shape (semicirc_num, 2)

    # Points along the diameter
    diam_x_points = np.linspace(-radius, radius, diam_num)
    diam_y_points = np.ones(diam_num) * (1 - radius)
    diam_points = np.array([diam_x_points, diam_y_points]).T  # Correct shaping (diam_num, 2)

    # Concatenate along axis 0 (rows)
    initial_points = np.concatenate((points, diam_points), axis=0)
    return initial_points


points_init = initial_points(0.975, 20, pathnum-20)
xinit = points_init[:,0]
yinit = points_init[:,1]

#Initial plot
main_tracks = []
acceleration_plot_tracks = []
x_tracks = []
y_tracks = []

fill = compute_fill(gam_0)


fixedPoint_plots = []
eos_plots = []
for i in range(0, pathnum):
    initial_conditions = [xinit[i], yinit[i]]  #Initial values for x and y
    

    lam = lambda_slide.get()
    gam = gamma_slide.get()

    # Solve the system of ODEs using odeint
    solution = odeint(ODEs, initial_conditions, N, args = (lam, gam))

    solution_x = solution[:, 0]
    solution_y = solution[:, 1]
    
    track_i = track_ax.plot(solution_x, solution_y, 'k', linewidth=.5, alpha=0.4)[0]
    main_tracks.append(track_i)
    track_i.set_visible(True)

    if i == 1:
        eos_plot_i, = gamma_ax.plot(N, gamma_phi(solution_x,solution_y), 'b-', linewidth=.5, alpha=0.4, label = r'$\gamma_\phi$')
    elif i % 1 == 0:
        eos_plot_i, = gamma_ax.plot(N, gamma_phi(solution_x,solution_y), 'b-', linewidth=.5, alpha=0.4)

    eos_plots.append(eos_plot_i)

    fixedPoints, fixedPoints_labels = fixedPoints_func(lam_0, gam_0)
    for point in fixedPoints:
        plot, = track_ax.plot(point[0], point[1], 'or', markeredgewidth=0.5, markeredgecolor='black')
        fixedPoint_plots.append(plot)


rScalingLine  = gamma_ax.plot([N[-1], N[0]], [4/3, 4/3], "r:", linewidth = .75, label='$\gamma_r=4/3$')
mScalingLine  = gamma_ax.plot([N[-1], N[0]], [1, 1], "g--", linewidth = .75, label='$\gamma_m=1$')
CCScalingLine = gamma_ax.plot([N[-1], N[0]], [0, 0], "k-.", linewidth = .75, label='$\gamma_\Lambda=0$')
gam = gamma_slide.get()
gam_ScaleLine, = gamma_ax.plot([N[-1], N[0]], [gam, gam], "m-", linewidth = .75, label=r'$\gamma$')


#_____________________Setting plot labels etc________________________________

track_ax.set_xlabel('$x$', x=1)
track_ax.set_ylabel('$y$', rotation = 0, y=1.01)
track_ax.set(xticks=[-1,-.5,0,.5,1], yticks=[0,.5,1],
             xlim = [-1.2,1.2], ylim=[-0.2,1.15],
             xticklabels = ['$-1$','$-1/2$','$0$','$1/2$','$1$'],
             yticklabels = ['$0$','$1/2$','$1$'])
track_ax.legend()

#Interactive plot, plot lines as in P40
lam_gam_ax.set_xlabel('$\lambda^2$', x=1.02)
lam_gam_ax.set_ylabel('$\gamma$', rotation = 0, y=1.02)
lam_gam_ax.set(xlim=[0,12], ylim=[0,2],
               yticks=[0,2/3,1,2], yticklabels=['$0$','$2/3$','$1$','$2$'],
               xticks=[0,2,6,12],  xticklabels=['$0$','$2$','$6$','$12$'])
lam_gam_ax.text(1,1.5,'Ia', size='14')
lam_gam_ax.text(3,1.5,'Ib', size='14')
lam_gam_ax.text(4,2/3,'II', size='14')
lam_gam_ax.text(9,1.5,'III', size='14')

lam_gam_ax.plot([0,6], [12,2],'k', linestyle = '-')
lam_gam_ax.plot([0,12],[0,0],'k')
lam_gam_ax.plot([0,0],[0,2],'k')
lam_gam_ax.plot([6,6],[0,2],'k')
lam_gam_ax.plot([0,6],[0,2],'k')
lam_gam_ax.plot([2,2],[0,2],'k', linestyle = ':')
#lam_gam_ax.plot([0, 12], [2/3, 2/3],'k', linestyle = ':')

fig.suptitle(f'$\lambda^2 = {lambda_slide.get()**2:.0f},\; \gamma={gamma_slide.get():.0f}$',fontsize=14)

gamma_ax.set(yticks = [0, 1, 4/3, 2], ylim=[-0.1,2.1],
            yticklabels = ['$0$','$1$', '$4/3$', '$2$'], xlim=[0,12])
gamma_ax.set_xlabel('$N$')
gamma_ax.set_ylabel('$\gamma_\phi$', rotation = 0)
gamma_ax.yaxis.set_ticks_position('both')
gamma_ax.tick_params(axis='x', which='both', labelbottom=True) 

#fig2.savefig("Figures/Arbitrary Fluid/EoS_lambda2_{}_gamma_{}.svg".format(round(lam_0**2),round(gam_0)), format='svg')
#fig3.savefig("Figures/Arbitrary Fluid/Track_lambda2_{}_gamma_{}.svg".format(round(lam_0**2),round(gam_0)), format='svg')

window.mainloop()
#End




