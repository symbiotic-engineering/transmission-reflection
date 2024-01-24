# testing hydrostatic calcs over different frequencies
import capytaine as cpt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.linalg import block_diag
from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity, froude_krylov_force

# initializing parameters
# wind turbines
rt = 4     # turbine radius [m]
lt = 50    # turbine length [m]
xt = 6000  # turbine 1 x-postion [m]
yt = 2000  # turbine 1 y-position [m]

# wave energy converters
rw = 10     # wec radius [m]
lw = 5      # wec length [m]
xw = 4000      # wec 1 x-position
yw = 2500      # wec 1 y-position
nr = 20         # number of panels along radius
ntheta = 20     # number of panels in theta direction
nz = 20         # number of panels in z-direction

# all bodies follow these parameters
z = 0      # all bodies z start position
w = 1.047        # frequency [rad/s]
B = 3*np.pi/4    # wave direction [rad]
depth = 40       # average water depth at southfork

# defining turbine mesh
turb = cpt.FloatingBody(mesh=cpt.mesh_vertical_cylinder(length=lt, radius=rt,center=(xt,yt,z),resolution=(nr,ntheta,nz),name='turb'))
turb.keep_immersed_part()
turb.center_of_mass=(0,0,-lt/2)
dofs = turb.add_all_rigid_body_dofs()
turb.inertia_matrix = turb.compute_rigid_body_inertia()
turb.hydrostatic_stiffness = turb.compute_hydrostatic_stiffness()

# define WEC mesh
wec = cpt.FloatingBody(mesh=cpt.mesh_vertical_cylinder(length=lw, radius=rw,center=(xw,yw,z),resolution=(nr,ntheta,nz),name='wec'))
wec.keep_immersed_part()
wec.center_of_mass=(0,0,-lw/2)
dofs = wec.add_all_rigid_body_dofs()
wec.inertia_matrix = wec.compute_rigid_body_inertia()
wec.hydrostatic_stiffness = wec.compute_hydrostatic_stiffness()

# create turbine array
t_array = turb + turb.translated((-2000,0,0),name='t2') + turb.translated((-4000,0,0),name='t3') + turb.translated((0,2400,0),name='t4') + turb.translated((-2000,2400,0),name='t5') + turb.translated((0,4800,0),name='t6') + turb.translated((-4000,4800,0),name='t7') + turb.translated((0,7200,0),name='t8') + turb.translated((-2000,7200,0),name='t9') + turb.translated((-4000,7200,0),name='t10') + turb.translated((-2000,9600,0),name='t11') + turb.translated((-4000,9600,0),name='t12')
t_array.add_all_rigid_body_dofs()
t_array.keep_only_dofs(dofs=['turb__Pitch','t2__Pitch','t3__Pitch','t4__Pitch','t5__Pitch','t6__Pitch','t7__Pitch','t8__Pitch','t9__Pitch','t10__Pitch','t11__Pitch','t12__Pitch'])

# create WEC array
w_array = wec + wec.translated((500,0,0),name='w2') + wec.translated((1000,0,0),name='w3') + wec.translated((1000,500,0),name='w4') + wec.translated((1000,1000,0),name='w5')
w_array.add_all_rigid_body_dofs()
w_array.keep_only_dofs(dofs=['wec__Heave','w2__Heave','w3__Heave','w4__Heave','w5__Heave'])

array = w_array + t_array

# solving hydrodynamics
solver = cpt.BEMSolver()
diff_prob = cpt.DiffractionProblem(body=array, wave_direction=B, water_depth=depth,omega=w)
diff_result = solver.solve(diff_prob,keep_details=(True))
rad_prob = [
    cpt.RadiationProblem(body=array, radiating_dof=dof, water_depth=depth,omega=1.047)
    for dof in array.dofs
    ]
rad_result = solver.solve_all(rad_prob,keep_details=(True))

dataset = cpt.assemble_dataset(rad_result + [diff_result])
#print('data',dataset)

# damping
damp = dataset['radiation_damping'].sel(radiating_dof=['wec__Heave','w2__Heave','w3__Heave','w4__Heave','w5__Heave'],
                                     influenced_dof=['wec__Heave','w2__Heave','w3__Heave','w4__Heave','w5__Heave'])
#print('damping',damp)
B = np.array([639909.55595804,640249.72070651,642644.96928585,639722.4077731,637023.25313503])
# added mass
add = dataset['added_mass'].sel(radiating_dof=['wec__Heave','w2__Heave','w3__Heave','w4__Heave','w5__Heave'],
                                     influenced_dof=['wec__Heave','w2__Heave','w3__Heave','w4__Heave','w5__Heave'])
#print('added mass',add)
A = np.array([1464351.62959333,1458179.71284097,1466053.99755036,1462288.77604565,1462356.96022831])
# hydrostatic stiffness
stiff = dataset['hydrostatic_stiffness'].sel(radiating_dof=['wec__Heave','w2__Heave','w3__Heave','w4__Heave','w5__Heave'],
                                     influenced_dof=['wec__Heave','w2__Heave','w3__Heave','w4__Heave','w5__Heave'])
#print('stiffness',stiff)
K = np.array([3031456.71481823,3031456.71481823,3031456.71481823,3031456.71481823,3031456.71481823])
# inertia matrix
mass = dataset['inertia_matrix'].sel(radiating_dof=['wec__Heave','w2__Heave','w3__Heave','w4__Heave','w5__Heave'],
                                     influenced_dof=['wec__Heave','w2__Heave','w3__Heave','w4__Heave','w5__Heave'])
#print('inertia matrix',mass)
M = np.array([772542.48593737,772542.48593737,772542.48593737,772542.48593737,772542.48593737])

# heave exciting forces
FK = [froude_krylov_force(diff_prob)['wec__Heave']]
FK += [froude_krylov_force(diff_prob)['w2__Heave']]
FK += [froude_krylov_force(diff_prob)['w3__Heave']]
FK += [froude_krylov_force(diff_prob)['w4__Heave']]
FK += [froude_krylov_force(diff_prob)['w5__Heave']]
#print('froude_krylov',FK)
FK = np.array(FK)

dif = [diff_result.forces['wec__Heave']]
dif += [diff_result.forces['w2__Heave']]
dif += [diff_result.forces['w3__Heave']]
dif += [diff_result.forces['w4__Heave']]
dif += [diff_result.forces['w5__Heave']]
#print('diffraction force',dif)
dif = np.array(dif)

ex_force = FK + dif
print('exciting force',ex_force)
F = np.array(ex_force)

omega = np.array([1.047, 1.047, 1.047, 1.047, 1.047])
#mag = np.abs(ex_force)
#print('total force',mag)
def wec_power(M, A, B, K, F, omega):
    # Calculates power absorbed by PTO in 2 different cases. 
    # M = mass of WEC; A = added mass; B = WEC damping coeff;
    # K = stiffness coeff; Bpto = pto damping coeff; Kpto = pto spring const
    # Fexc = wave excitation force; omega = incoming wave frequency;

    
    # Case 1: Kpto & Bpto (reactive control)
    
    # Define optimal Bpto and Kpto
    Bpto_1 = B
    Kpto_1 = omega**2*(M+A)-K  
    
      
    # WEC motion (complex) 
    X_1 = F/((-omega**2)*(M+A)+(B+Bpto_1)*omega*1j+K+Kpto_1)    
    # WEC power case 1
    P_1 = 0.5*Bpto_1*(abs(X_1))**2*omega**2

    # Case 2: Bpto only (passive control)
    
    # For Kpto = 0, define optimal Bpto
    Kpto_2 = 0
    Bpto_2 = (B**2+(omega*(M+A)-K/omega)**2)**0.5
    
    # wec motion
    X_2 = F/((-omega**2)*(M+A)+(B+Bpto_2)*omega*1j+K+Kpto_2)
    # power case 2
    P_2 = 0.5*Bpto_1*(abs(X_2))**2*omega**2
         
    return[P_1, P_2]

[p1, p2]=wec_power(M, A, B, K, F, omega)

print('power 1',p1)
print('power 2',p2)

#def plotting(diff_result,xw,yw,xt,yt):
# plotting the free surface pertubation 

# generate the free surface
free_surface = cpt.FreeSurface(x_range=(-2000, 10000), y_range=(-2000, 10000), nx=200, ny=200)
dif_el = solver.get_free_surface_elevation(diff_result, free_surface)
# add incoming waves

h_i = free_surface.incoming_waves(diff_result)    # not double counting diffraction elevation from spheres. diff_result stores airys wave data
kt = dif_el/h_i
#h_t = (dif_el + h_i)                              # total perturbed wave field
#kd = h_t/h_i                                      # disturbance coefficient

# plots
x = np.linspace(-2000,10000,200)
y = np.linspace(-2000,10000,200)
X, Y = np.meshgrid(x, y)
Z = kt.reshape(200,200)
fig, ax = plt.subplots()
CS = ax.contourf(X,Y,Z,cmap="viridis")
plt.plot(yw,xw,'ro',markersize=6)
plt.plot(yw,xw+500,'ro',markersize=6)
plt.plot(yw,xw+1000,'ro',markersize=6)
plt.plot(yw+500,xw+1000,'ro',markersize=6)
plt.plot(yw+1000,xw+1000,'ro',markersize=6)

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
#return[h_i,dif_el]

plt.show()