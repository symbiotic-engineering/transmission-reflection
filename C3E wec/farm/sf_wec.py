# testing hydrostatic calcs over different frequencies
import capytaine as cpt
from scipy.linalg import block_diag
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.linalg import block_diag
from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity, froude_krylov_force

# initializing parameters
r = 10     # radius [m]
l = 5
x = 0
y = 0
nr = 20         # number of panels along radius
ntheta = 20     # number of panels in theta direction
nz = 20         # number of panels in z-direction
z = 0
w = 1.047        # frequency [rad/s]
B = 3*np.pi/4      # wave direction [rad]
depth = 40      # average water depth at southfork

# defining mesh
body = cpt.FloatingBody(mesh=cpt.mesh_vertical_cylinder(length=l, radius=r,center=(x,y,z),resolution=(nr,ntheta,nz),name='cyl'))
body.keep_immersed_part()
body.center_of_mass=(0,0,-l/2)
dofs = body.add_all_rigid_body_dofs()
body.inertia_matrix = body.compute_rigid_body_inertia()
body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()

# create array
array = body + body.translated((40,40,0),name='2') + body.translated((75,75,0),name='3') + body.translated((150,150,0),name='4') + body.translated((200,200,0),name='5')
from scipy.linalg import block_diag
array.add_all_rigid_body_dofs()
array.keep_only_dofs(dofs=['cyl__Heave','2__Heave','3__Heave','4__Heave','5__Heave'])
array.show_matplotlib(normal_vectors = True)

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

# damping
damp = dataset['radiation_damping'].sel(radiating_dof=['cyl__Heave','2__Heave','3__Heave','4__Heave','5__Heave'],
                                     influenced_dof=['cyl__Heave','2__Heave','3__Heave','4__Heave','5__Heave'])
print('damping',damp)
B = np.array([631508.66601067,632501.09970108,629862.74410561,632501.09970105,631508.66601068])
# added mass
add = dataset['added_mass'].sel(radiating_dof=['cyl__Heave','2__Heave','3__Heave','4__Heave','5__Heave'],
                                     influenced_dof=['cyl__Heave','2__Heave','3__Heave','4__Heave','5__Heave'])
print('added mass',add)
A = np.array([1465611.97398437,1472291.51181506,1469126.68547142,1472291.51181503,1465611.97398438])
# hydrostatic stiffness
stiff = dataset['hydrostatic_stiffness'].sel(radiating_dof=['cyl__Heave','2__Heave','3__Heave','4__Heave','5__Heave'],
                                     influenced_dof=['cyl__Heave','2__Heave','3__Heave','4__Heave','5__Heave'])
print('stiffness',stiff)
K = np.array([3031456.71481823,3031456.71481823,3031456.71481823,3031456.71481823,3031456.71481823])
# inertia matrix
mass = dataset['inertia_matrix'].sel(radiating_dof=['cyl__Heave','2__Heave','3__Heave','4__Heave','5__Heave'],
                                     influenced_dof=['cyl__Heave','2__Heave','3__Heave','4__Heave','5__Heave'])
print('inertia matrix',mass)
M = np.array([772542.48593737,772542.48593737,772542.48593737,772542.48593737,772542.48593737])

# heave exciting forces
FK = [froude_krylov_force(diff_prob)['cyl__Heave']]
FK += [froude_krylov_force(diff_prob)['2__Heave']]
FK += [froude_krylov_force(diff_prob)['3__Heave']]
FK += [froude_krylov_force(diff_prob)['4__Heave']]
FK += [froude_krylov_force(diff_prob)['5__Heave']]
#print('froude_krylov',FK)
FK = np.array(FK)

dif = [diff_result.forces['cyl__Heave']]
dif += [diff_result.forces['2__Heave']]
dif += [diff_result.forces['3__Heave']]
dif += [diff_result.forces['4__Heave']]
dif += [diff_result.forces['5__Heave']]
#print('diffraction force',dif)
dif = np.array(dif)

ex_force = FK + dif
print('exciting force',ex_force)
F = np.array(ex_force)

# omega = np.array([1.047, 1.047, 1.047, 1.047, 1.047])
# #mag = np.abs(ex_force)
# #print('total force',mag)
# def wec_power(M, A, B, K, F, omega):
#     # Calculates power absorbed by PTO in 2 different cases. 
#     # M = mass of WEC; A = added mass; B = WEC damping coeff;
#     # K = stiffness coeff; Bpto = pto damping coeff; Kpto = pto spring const
#     # Fexc = wave excitation force; omega = incoming wave frequency;

    
#     # Case 1: Kpto & Bpto (reactive control)
    
#     # Define optimal Bpto and Kpto
#     Bpto_1 = B
#     Kpto_1 = omega**2*(M+A)-K  
    
      
#     # WEC motion (complex) 
#     X_1 = F/((-omega**2)*(M+A)+(B+Bpto_1)*omega*1j+K+Kpto_1)    
#     # WEC power case 1
#     P_1 = 0.5*Bpto_1*(abs(X_1))**2*omega**2

#     # Case 2: Bpto only (passive control)
    
#     # For Kpto = 0, define optimal Bpto
#     Kpto_2 = 0
#     Bpto_2 = (B**2+(omega*(M+A)-K/omega)**2)**0.5
    
#     # wec motion
#     X_2 = F/((-omega**2)*(M+A)+(B+Bpto_2)*omega*1j+K+Kpto_2)
#     # power case 2
#     P_2 = 0.5*Bpto_2*(abs(X_2))**2*omega**2
         
#     return[P_1, P_2]

# [p1, p2]=wec_power(M, A, B, K, F, omega)

# print('power 1',p1)
# print('power 2',p2)

# # post-processing
# # creating mesh of free surface
# free_surface = cpt.FreeSurface(x_range=(-500, 1500), y_range=(-500, 1500), nx=400, ny=400)

# dif_el = solver.get_free_surface_elevation(diff_result, free_surface)
# #rad_el = solver.get_free_surface_elevation(rad_result, free_surface)
# # add incoming waves

# h_i = free_surface.incoming_waves(diff_result)    # not double counting diffraction elevation from spheres. diff_result stores airys wave data
# h_t = (dif_el + h_i)                              # total perturbed wave field
# for i in range(0, len(free_surface.mesh.faces_centers)):
#     free_surface.mesh.faces_centers[i][2] = h_t[i]

# kd = h_t/h_i             # disturbance coefficient

# kt = dif_el/h_i          # transmission coeff

# # plots
# x = np.linspace(-500,1500,400)
# y = np.linspace(-500,1500,400)
# X, Y = np.meshgrid(x, y)
# Z = kt.reshape(400,400)
# fig, ax = plt.subplots()
# CS = ax.contourf(X,Y,Z,cmap="viridis",vmin = 0, vmax = 1)
# plt.plot(0,0,'ro',markersize=6)
# plt.plot(250,250,'ro',markersize=6)
# plt.plot(500,500,'ro',markersize=6)
# plt.plot(750,750,'ro',markersize=6)
# plt.plot(1000,1000,'ro',markersize=6)

# formatter = ticker.ScalarFormatter(useMathText=True)
# formatter.set_scientific(True)
# formatter.set_powerlimits((-1,1))
# ax.yaxis.set_major_formatter(formatter)
# ax.xaxis.set_major_formatter(formatter)

# cbar = fig.colorbar(CS)
# cbar.set_label('$K_{t/r}$',fontsize=17,rotation=270,labelpad=17)
# ax.set_xlabel('x [m]',fontsize=17)
# ax.set_ylabel('y [m]',fontsize=17)
# cbar.ax.tick_params(labelsize=15)
# plt.xticks(fontsize=14, rotation=90)
# plt.yticks(fontsize=14, rotation=90)
# fig.tight_layout()
# plt.savefig('rm3array_wave_transimission.pdf')
# plt.show()