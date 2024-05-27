#_______________________________Imports______________________________________
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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




#_________________________Set up figures________________________________
fig2 = plt.figure(figsize=(6.5, 4.5))
fig2.set_facecolor('white')

fig3 = plt.figure(figsize=(9, 6))
fig3.set_facecolor('white')

fig5 = plt.figure(figsize=(7.5, 5))
fig5.set_facecolor('white')

fig6 = plt.figure(figsize=(9, 6))
fig6.set_facecolor('white')

fig7 = plt.figure(figsize=(7, 4.5))
fig7.set_facecolor('white')



track_axis_dims2 = [0,.075,.9,.9]
track_ax = fig2.add_axes(track_axis_dims2, projection='3d')
track_ax.view_init(elev=24, azim=66)

cbar_ax_dims2 = [.7,.25,.02,.6]
cbar_ax = fig2.add_axes(cbar_ax_dims2)

dens_axis_dims2 = [.1,.15,.825,.75]
dens_ax = fig3.add_axes(dens_axis_dims2)

gamma_axis_dims2 = [.1,.15,.8,.75]
gamma_ax = fig5.add_axes(gamma_axis_dims2)

hubble_ax_dims2 = [.15,.15,.8,.8]
hubble_ax = fig6.add_axes(hubble_ax_dims2)

rho_ax_dims = [.15,.15,.8,.8]
rho_ax = fig7.add_axes(rho_ax_dims)

figures = [fig2, fig3, fig5, fig6, fig7]


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

lam_0s = [[0,0], [1,-1], [10, 0.1], [30,-30]]#, [20, 0.2]]

lam1_0 = -0.4
lam1_min = -30
lam1_max = 30

lam2_0 = 4
lam2_min = -30
lam2_max = 30

lam1_0, lam2_0 = lam_0s[3]


Ni = 25

N = np.linspace(0, Ni, abs(int(Ni*10000)))

c = 1 #3e5     # Given in km/s
##V = z * c   # ""
h = 0.738
H_0 = 100*h # km/(s Mpc)
xinit = np.linspace(-0.99, 0.99, pathnum)

H0_SN = 73.04
H0_SN_Err = 1.04

H0_PL = 68.75
H0_PL_Err = 0.52

#Observed Values
Omega_phi_0 = 0.738
Omega_Lambda0 = 0.738
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

#For lam1,lam2 = 1,-1 use Nplotmax = 4, x=2e-8, y1,y2 = 1e-10, z=0.995
#For lam1=10, lam2 = 0.1, use Nplotmax = 4, x=2e-5, y1=1e-5 y2=1e-10 z=0.99
N_plot_max = 5
N_plot_min = -15

# state_0 = [0.163,
#           0.115,
#           0.0000000001,
#           0.979]

##state_0 = [0.00002,
##           0.00001,
##           0.0000000001,
##           0.99]

state_0 = [0.16329,
           0.11547,
           2.6735e-15,
           0.97979]

state_0 = [0,
           1e-12,
           1e-12,
           0.999796]

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

        if lambdaValue >= 0: #Makes sure to only plot positive y points
            # C1 and C2
            if lambdaValue**2 > 3: #Makes sure to only plot positive y points
                fixedPoints = np.append(fixedPoints,
                                        [[np.sqrt(3/2)/lambdaValue,
                                            np.sqrt(3/2)/lambdaValue, 0]], axis=0)
                fixedPoints_labels.append('$C' + str(n+1) + '$')

            # E1 and E2
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

#___________________________________Initial Plot_____________________________

#Initial plot
fixedPoint_plots = []



lam1 = lam1_0
lam2 = lam2_0

# Solve the system of ODEs using odeint
solution = odeint(ODEs, state_0, N, args = (np.array([lam1, lam2]),))
pathx  = solution[:, 0]
pathy1 = solution[:, 1]
pathy2 = solution[:, 2]
pathz  = solution[:, 3]

pathy = np.sqrt(pathy1**2 + pathy2**2)

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
m_max_val = getRedshift(NAxis[indexMPeak])
msf_eq_val = getRedshift(NAxis[indexMPhi_eq])


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

rScalingLine,  = gamma_ax.plot([N[-1], N[0]], [4/3, 4/3], "r:", linewidth = .75, label='$\gamma_r=4/3$')
mScalingLine,  = gamma_ax.plot([N[-1], N[0]], [1, 1], "g--", linewidth = .75, label='$\gamma_m=1$')
accel_plot, = gamma_ax.plot(NAxis,
                accelerationExpression(pathx,pathy,pathz),'darkorange', label='Acceleration')
effective_eos, = gamma_ax.plot(NAxis, gamma_phi(pathx, pathy), 'b-', linewidth=1, label = r'$\gamma_\phi$')
CCScalingLine, = gamma_ax.plot([N[-1], N[0]], [0, 0], "k-.", linewidth = .75, label='$\gamma_\Lambda=0$')

#bcgd_gam = (mass_dens + rad_dens*4/3)/2
#backgrd_scaling, = gamma_ax.plot(NAxis, bcgd_gam, "m-", linewidth = .75, label=r'$(\gamma_r\Omega_r+\gamma_m\Omega_m)/2$')

# y1_dens_plot, =  dens_ax.plot(NAxis, pathy1**2, 'b--',
#        label = "$y_1^2$")
# y2_dens_plot, =  dens_ax.plot(NAxis, pathy2**2, 'b--',
#        label = "$y_2^2$")



x_i, y_i, z_i = state_0[0], np.sqrt(state_0[1]**2 + state_0[2]**2), state_0[3]


# state0_point, = track_ax.plot(x_i,y_i,z_i, 'cX')
track = track_ax.plot(
                pathx, pathy, pathz, 'b', linewidth=2)[0]
fixedPoints, fixedPoints_labels = fixedPoints_func([lam1,lam2])
for point in fixedPoints:
    plot, = track_ax.plot(point[0], point[1], point[2], 'or', markeredgewidth=0.5, markeredgecolor='black')
    fixedPoint_plots.append(plot)


track.set_visible(True)


xRunningIntegral = np.append(0, cumulative_trapezoid(pathx, NAxis))

hubbleFromY1 = HfromY(pathy1, xRunningIntegral, lam1)
hubbleFromY1 *= H0_PL / hubbleFromY1[indexToday]

hubbleFromY2 = HfromY(pathy2, xRunningIntegral, lam2)
hubbleFromY2 *= H0_PL / hubbleFromY2[indexToday]

rho_r = rad_dens * hubbleFromY2**2
rho_m = mass_dens * hubbleFromY2**2
rho_phi = phi_dens * hubbleFromY2**2


z_decouple = 1100
N_decouple = -np.log(1+z_decouple)
rho_ax.fill_betweenx([1e-0,1e20], N_plot_min, NAxis[indexMR_eq], alpha=0.1, color="pink")
rho_ax.fill_betweenx([1e-0,1e20], NAxis[indexMR_eq], NAxis[indexMPhi_eq], alpha=0.1, color="lime")
rho_ax.fill_betweenx([1e-0,1e20], NAxis[indexMPhi_eq], N_plot_max, alpha=0.1, color="cyan")

rho_ax.plot([NAxis[indexMR_eq],NAxis[indexMR_eq]],[1e-0,1e20],'m:', linewidth=1.2, label=r'$\rho_r=\rho_m$')
rho_ax.plot([NAxis[indexMPhi_eq],NAxis[indexMPhi_eq]],[1e-0,1e20],'c-.', linewidth=1.2, label=r'$\rho_m=\rho_\phi$')
rho_ax.plot([N_decouple,N_decouple],[1e0,1e20],'darkorange', linestyle='--', linewidth=1.2, label=r'$z_{\mathrm{dec}}\approx 1100$')
rho_r_plot, = rho_ax.plot(NAxis, rho_r, "r", label = r"$\rho_r\sim a^{-4}$")
rho_m_plot, = rho_ax.plot(NAxis, rho_m, "g", label = r"$\rho_m\sim a^{-3}$")
rho_phi_plot, = rho_ax.plot(NAxis, rho_phi, "b", label = r"$\rho_\phi \sim 1$")



N_start_idx = np.argmax(N > N_plot_min)
N_end_idx   = np.argmax(N > N_plot_max)


hubble_ax.plot([0,0],[1e0,1e20],'k-', linewidth=1, label='Today')
hubble_ax.plot([N_decouple,N_decouple],[1e0,1e20],'darkorange', linestyle=':', linewidth=1, label=r'$z_{\mathrm{dec}}\approx 1100$')
hubble_plot, = hubble_ax.plot(NAxis[N_start_idx:N_end_idx], hubbleFromY2[N_start_idx:N_end_idx], color='b',label=r'$H^{(\phi)}(N)$')
#hubble_plot, = hubble_ax.plot(NAxis, hubbleFromY1)

hubble_SN_line = hubble_ax.plot([N_plot_min,N_plot_max], [H0_SN, H0_SN], alpha=0.75, color = "cyan", label=r'$H_0^{\mathrm{SN}}=73.04\pm 1.04$')
hubble_PL_line = hubble_ax.plot([N_plot_min,N_plot_max], [H0_PL, H0_PL], alpha=0.75, color = "magenta", label=r'$H_0^{\mathrm{PL}}=67.85\pm 0.52$')






'''Figure Production Settings'''
track_ax.set(xlabel='$x$', ylabel='$y$', zlabel='$z$',
             xlim = [-1,1], ylim = [0,1], zlim = [0,1],
             xticks = [-1, -0.5, 0, 0.5, 1],
             yticks = [0, 0.5, 1],
             zticks = [0, 0.5, 1])
track_ax.set_box_aspect([2, 1, 1])
track_ax.axis("off")



# static_line = accel_ax.plot([N[-1], N[0]],[0,0], "k--", linewidth = 0.5)
# accel_ax.set(ylabel="Acceleration", ylim=[-1.1,1.1],
#              yticks=[-1,-1/2,0,1/2,1], yticklabels = ['$-1$','$-1/2$', '$0$', '$1/2$', '$1$'],
#                xlim=[-8,3], xticks = [-8,-6,-4,-2,0,2], xlabel="$N$",
#             xticklabels = ['$-8$', '$-6$', '$-4$', '$-2$','$0$','$2$'])
# accel_ax.yaxis.set_ticks_position('both')
# accel_ax.tick_params(axis='x', which='both', labelbottom=True)


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


dens_ax.text(NAxis[indexMR_eq]-1,.1,f'$z={mr_eq_val:.1f}$',backgroundcolor='1', fontsize=12)
dens_ax.text(NAxis[indexMPeak]-1.2,1.1,f'$z={m_max_val:.1f}$',backgroundcolor='1', fontsize=12)
dens_ax.text(NAxis[indexMPhi_eq]-1.3,1.1,f'$z={msf_eq_val:.1f}$',backgroundcolor='1', fontsize=12)

dens_ax.set(xlabel="$N$", ylabel="Density Parameters",
            ylim=[-0.1,1.2],yticks=[0,1/4,1/2,3/4,1],
            yticklabels = ['$0$','$1/4$','$1/2$','$3/4$','$1$'],
            xlim = [N_plot_min, N_plot_max], xticks=[-15,-10,-5,0,5])




# todayLine, = dens_ax.plot([0,0], [-0.2,1.2], 'k', label='Today:\; $z = 0$')
# MR_eqLine, = dens_ax.plot([NAxis[indexMR_eq],NAxis[indexMR_eq]], [-0.2,1.2], 'k:', linewidth = 0.75,
#                           label = f'$\Omega_m=\Omega_r:\; z={mr_eq_val:.3f}$')
# MPeakLine, = dens_ax.plot([NAxis[indexMPeak],NAxis[indexMPeak]], [-0.2,1.2], 'k--', linewidth = 0.75,
#                           label = f'max$(\Omega_m):\; z={m_max_val:.3f}$')
# MPhi_eqLine, = dens_ax.plot([NAxis[indexMPhi_eq],NAxis[indexMPhi_eq]], [-0.2,1.2], 'k-.', linewidth = 0.75,
#                           label = f'$\Omega_m=\Omega_\phi:\; z={msf_eq_val:.3f}$')

rho_ax.set_yscale('log', base=10, subs=[10**x
                         for x in (0.25, 0.5, 0.75)], nonpositive='mask')
rho_ax.set( xlabel = "$N$",
            ylabel=r"$\log_{10}(\rho)$",
            xlim = [N_plot_min, N_plot_max],
            ylim=[1e0,1e20],
            yticks=[1e0,1e5,1e10,1e15,1e20],
            yticklabels=['$0$','$5$','$10$', '$15$',
                             '$20$'],
            xticks=[-15,-10,-5,0,5])
rho_ax.legend(loc= 'upper right', ncol=2) 





hubble_ax.set_yscale('log', base=10, subs=[10**x
                         for x in (0.25, 0.5, 0.75)], nonpositive='mask')
hubble_ax.set(  xlabel = "$N$",
                ylabel=r"$\log_{10}(H)$",
                xlim = [N_plot_min, N_plot_max],
                ylim = [1e0,1e15],
                yticks=[1e0,1e5,1e10,1e15],
                yticklabels=['$0$', '$5$', '$10$','$15$'],
                xticks=[-15,-10,-5,0,5])
hubble_ax.legend(loc='upper left')

axins = hubble_ax.inset_axes([0.5,0.5,0.45,0.45])
axins.set_xlim(-1, 2)
axins.set_ylim(50,80)
hubble_ax.indicate_inset_zoom(axins, edgecolor="black")
axins.set(xticks=[-1,0,1,2], yticks=[50,60,70,80])
axins.tick_params(axis='both', which='major', labelsize=9)


N_start_idx_inset = np.argmax(N > -1)
N_end_idx_inset   = np.argmax(N > 2)

hubble_ax.plot([0,0],[50,80],'k-', linewidth=1)
axins.plot(NAxis[N_start_idx_inset:N_end_idx_inset], 
           hubbleFromY2[N_start_idx_inset:N_end_idx_inset],'b')
axins.plot([N_plot_min,N_plot_max], [H0_SN, H0_SN], alpha=0.75, color = "cyan")
axins.plot([N_plot_min,N_plot_max], [H0_PL, H0_PL], alpha=0.75, color = "magenta")
axins.fill_between(NAxis[N_start_idx_inset:N_end_idx_inset], H0_SN + H0_SN_Err, H0_SN - H0_SN_Err, alpha=0.2, color="cyan")
axins.fill_between(NAxis[N_start_idx_inset:N_end_idx_inset], H0_PL + H0_PL_Err, H0_PL - H0_PL_Err, alpha=0.2, color="magenta")



#Additional code for making paper plots
legend_lines1 = []
legend_lines1.append([todayLine, MR_eqLine, MPeakLine, MPhi_eqLine])
legend_lines2 = []
legend_lines2.append([Radn_dens_plot, Mass_dens_plot, Phi_dens_plot])

legend_lines_rho = []
legend_lines_rho.append([rho_r_plot, rho_m_plot, rho_phi_plot])


legend1 = dens_ax.legend(legend_lines1[0], ["Today","$\Omega_m=\Omega_r$","max$(\Omega_m)$",
                                            "$\Omega_m=\Omega_\phi$"], loc='upper left', fontsize=12)
legend2 = dens_ax.legend(legend_lines2[0], ['$\Omega_r$', '$\Omega_m$', '$\Omega_\phi$'],
                         loc='center left', bbox_to_anchor=(0, .45), fontsize=12)
dens_ax.add_artist(legend1)
dens_ax.yaxis.set_ticks_position('both')




##d_lum_ax.set(ylabel = "$d_L$", xlabel= '$z$',
##              xlim=[0,3], ylim=[0,6],
##              xticks=[0,1,2,3], yticks=[0,1,2,3,4,5,6],
##              xticklabels = ['$0$','$1$','$2$', '$3$'],
##              yticklabels = ['$0$','$1$','$2$', '$3$', '$4$','$5$','$6$'])
##d_lum_ax.legend(loc=4)



#fig2.savefig("Figures/Two Potentials/UnlabeledTrack_Near_B.svg", format='svg')
#fig3.savefig("Figures/Two Potentials/Density_LCDM_Approximation.svg", format='svg')
#fig4.savefig("Figures/Two Potentials/Accel_lambda1_{}_lambda2_{}.svg".format(round(lam1_0),round(lam2_0)), format='svg')
#fig5.savefig("Figures/Two Potentials/Gamma_Near_B.svg", format='svg')
#fig6.savefig("Figures/Two Potentials/Rho_LCDM_Approximation.svg", format='svg') 
#fig6.savefig("Figures/Two Potentials/Hubble_Near_B.svg", format='svg') 

plt.show()




