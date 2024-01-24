# testing hydrostatic calcs over different frequencies
import capytaine as cpt
from capytaine.bem.airy_waves import froude_krylov_force
from scipy.linalg import block_diag
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

# initializing parameters
r = 4     # radius [m]
l = 50
x = 6000
y = 2000
z = 0
w = 1.047        # frequency [rad/s]
B = 3*np.pi/4    # wave direction [rad]
depth = 40       # average water depth at southfork

# defining mesh
body = cpt.VerticalCylinder(length=l, radius=r,center=(x,y,z),name='cyl')
#body = cpt.bodies.predefined.spheres.Sphere(radius=r,center=(x,y,z),ntheta=40,nphi=40,clip_free_surface=True)
body.keep_immersed_part()
body.center_of_mass=(0,0,-l/2)
dofs = body.add_all_rigid_body_dofs()
body.inertia_matrix = body.compute_rigid_body_inertia()
body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()

# create array
array = body + body.translated((-2000,0,0),name='2') + body.translated((-4000,0,0),name='3') + body.translated((0,2400,0),name='4') + body.translated((-2000,2400,0),name='5') + body.translated((0,4800,0),name='6') + body.translated((-4000,4800,0),name='7') + body.translated((0,7200,0),name='8') + body.translated((-2000,7200,0),name='9') + body.translated((-4000,7200,0),name='10') + body.translated((-2000,9600,0),name='11') + body.translated((-4000,9600,0),name='12')
array.add_all_rigid_body_dofs()
array.keep_only_dofs(dofs=['cyl__Pitch','2__Pitch','3__Pitch','4__Pitch','5__Pitch','6__Pitch','7__Pitch','8__Pitch','9__Pitch','10__Pitch','11__Pitch','12__Pitch'])

# solving hydrodynamics
solver = cpt.BEMSolver()
diff_prob = cpt.DiffractionProblem(body=array, wave_direction=B, water_depth=depth, omega=w)
diff_result = solver.solve(diff_prob,keep_details=(True))
rad_prob = [
    cpt.RadiationProblem(body=array, radiating_dof=dof, water_depth=depth,omega=1.047)
    for dof in array.dofs
    ]
rad_result = solver.solve_all(rad_prob,keep_details=(True))

# post-processing
# creating mesh of free surface
free_surface = cpt.FreeSurface(x_range=(0, 12000), y_range=(0, 12000), nx=100, ny=100)

dif_el = solver.get_free_surface_elevation(diff_result, free_surface)
#rad_el = solver.get_free_surface_elevation(rad_result, free_surface)
# add incoming waves

h_i = free_surface.incoming_waves(diff_result)    # not double counting diffraction elevation from spheres. diff_result stores airys wave data
h_t = (dif_el + h_i)     # total perturbed wave field
kd = h_t/h_i                                            # disturbance coefficient

# plots
x = np.linspace(0,12000,100)
y = np.linspace(0,12000,100)
X, Y = np.meshgrid(x, y)
Z = kd.reshape(100,100)
fig, ax = plt.subplots()
CS = ax.contourf(X,Y,Z,cmap="viridis")
plt.plot(2000,6000,'ro',markersize=6)
plt.plot(2000,4000,'ro',markersize=6)
plt.plot(2000,2000,'ro',markersize=6)
plt.plot(4400,6000,'ro',markersize=6)
plt.plot(4400,4000,'ro',markersize=6)
plt.plot(6800,6000,'ro',markersize=6)
plt.plot(6800,2000,'ro',markersize=6)
plt.plot(9200,6000,'ro',markersize=6)
plt.plot(9200,4000,'ro',markersize=6)
plt.plot(9200,2000,'ro',markersize=6)
plt.plot(11600,4000,'ro',markersize=6)
plt.plot(11600,2000,'ro',markersize=6)

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-1,1))
ax.yaxis.set_major_formatter(formatter)
ax.xaxis.set_major_formatter(formatter)

cbar = fig.colorbar(CS)
cbar.set_label('$K_{D}$',fontsize=17,rotation=270,labelpad=17)
ax.set_xlabel('x [m]',fontsize=17)
ax.set_ylabel('y [m]',fontsize=17)
cbar.ax.tick_params(labelsize=15)
plt.xticks(fontsize=14, rotation=90)
plt.yticks(fontsize=14, rotation=90)
#ax.annotate(r'$\beta$ = $\pi$/2',xy=(-100,175),xytext=(-185,175),arrowprops=dict(arrowstyle='->'))
#ax.annotate(r'$\omega$ = 1.047 rad/s',xy=(-185,145))
plt.show()
#plt.savefig('T=6s.pdf')