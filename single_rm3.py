# testing hydrostatic calcs over different frequencies
import capytaine as cpt
from scipy.linalg import block_diag
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.linalg import block_diag
from scipy.interpolate import griddata
from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity, froude_krylov_force
from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
import pandas as pd
import csv

# initializing parameters
r = 10     # radius [m]
l = 5
x = 0
y = 0
nr = 40         # number of panels along radius
ntheta = 15     # number of panels in theta direction
nz = 10         # number of panels in z-direction
z = 0
w = 0.50       # frequency [rad/s]
B = np.pi/2      # wave direction [rad]
depth = 40      # average water depth at southfork
def wave_num(w):
    # for infinite depth
    g = 9.81            # gravitational constant (m/s^2)
    k_inf = w**2/g      # wave number infinite depth (rad^2/m)
    lam_inf = 2*np.pi/k_inf    # wavelength infinite depth (m)
    return k_inf, lam_inf

k, lambda_wave = wave_num(w)

# defining mesh
body = cpt.FloatingBody(mesh=cpt.mesh_vertical_cylinder(length=l, radius=r,center=(x,y,-1/2),resolution=(nr,ntheta,nz),name='cyl'))
body.keep_immersed_part()
body.center_of_mass=(0,0,-l/2)   #-l/2)
body.add_all_rigid_body_dofs()
body.inertia_matrix = body.compute_rigid_body_inertia()
body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
body.keep_only_dofs(dofs='Heave')
# body.show_matplotlib(normal_vectors = True)
# solving hydrodynamics
solver = cpt.BEMSolver()
diff_prob = cpt.DiffractionProblem(body=body, wave_direction=B, water_depth=depth,omega=w)
diff_result = solver.solve(diff_prob,keep_details=(True))

rad_prob = [
   cpt.RadiationProblem(body=body, radiating_dof=dof, water_depth=depth,omega=w)
   for dof in body.dofs
   ]
rad_result = solver.solve_all(rad_prob,keep_details=(True))

dataset = cpt.assemble_dataset(rad_result + [diff_result])

# hydro coeffs
damp = dataset['radiation_damping'].sel(radiating_dof=['Heave'],
                                     influenced_dof=['Heave'])
add = dataset['added_mass'].sel(radiating_dof=['Heave'],
                                     influenced_dof=['Heave'])
#stiff = dataset['hydrostatic_stiffness'].sel(radiating_dof=['Heave'],
#                                     influenced_dof=['Heave'])
#stiffnes_matrix = body.hydrostatic_stiffness.values
#mass = dataset['inertia_matrix'].sel(radiating_dof=['Heave'],
#                                     influenced_dof=['Heave'])
# print('added mass',add)
# print('damping',damp)
# print('stiffness',stiff)
# print('mass',mass)

# from capytaine cookbook
x1 = int(-2*lambda_wave)
x2 = int(2*lambda_wave)
nx = int(lambda_wave)
y1 = int(-2*lambda_wave)
y2 = int(2*lambda_wave)
ny = int(lambda_wave)
grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
diffraction = solver.compute_free_surface_elevation(grid, diff_result)
chosen_rad_result = rad_result[0]
radiation = solver.compute_free_surface_elevation(grid, chosen_rad_result)
incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)
total = diffraction + radiation + incoming_fse
kd = total/incoming_fse
np.savetxt('kd_050.csv', kd, delimiter=',')
np.savetxt('dif_plus_inc_data_bw.csv', total, delimiter=",")
np.savetxt('incoming_fse_data_bw.csv',incoming_fse,delimiter=',')

Z = np.real(total)
X = grid[0]
Y = grid[1]
plt.pcolormesh(X, Y, Z, cmap='inferno')
plt.xlabel("x")
plt.ylabel("y")
colorbar = plt.colorbar()
colorbar.set_label('Elevation')
#colorbar.set_ticks([-1, 0, 1])  # Set custom ticks
plt.tight_layout()
plt.savefig(f'elev_disturb/050.pdf')
plt.show()
#np.savetxt('elevation_data_rm3.csv', Z, delimiter=",")