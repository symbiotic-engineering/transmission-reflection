# testing hydrostatic calcs over different frequencies
import capytaine as cpt
from scipy.linalg import block_diag
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.linalg import block_diag
from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity, froude_krylov_force
from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
import pandas as pd
import csv

# initializing parameters
wi = 25         # width of flap [m]
th = 1          # thickness of flap [m]
h = 19          # height of flap [m], draft = 16 m
x = 0
y = 0
draft = 16 # m
z = 0.5*h-draft # box center [m]
cog = -11.4     # center of gravity [m]
nw = 20         # number of panels along width (x)
nt = 5          # number of panels along thickness (y)
nh = 20         # number of panels along height (z)
w = 0.50        # frequency [rad/s]
B = 0    #3*np.pi/4      # wave direction [rad]
depth = 40      # average water depth at southfork
def wave_num(w):
    # for infinite depth
    g = 9.81            # gravitational constant (m/s^2)
    k_inf = w**2/g      # wave number infinite depth (rad^2/m)
    lam_inf = 2*np.pi/k_inf    # wavelength infinite depth (m)
    return k_inf, lam_inf

k, lambda_wave = wave_num(w)

# defining mesh
body = cpt.FloatingBody(cpt.meshes.predefined.rectangles.mesh_parallelepiped(size=(th, wi, h), 
                                                                             resolution=(nt, nw, nh), 
                                                                             center=(x, y, z),
                                                                             name='rect'))
body.keep_immersed_part()
body.center_of_mass=(0,0,cog)
dofs = body.add_all_rigid_body_dofs()
body.inertia_matrix = body.compute_rigid_body_inertia()
body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
body.keep_only_dofs(dofs='Pitch')
# # Check:
# awl = wi * th  # m**2
# kb = 0.5 * draft  # m
# kg = draft+cog # m
# bmt = wi ** 2 / 12 / draft  # m
# bml = th ** 2 / 12 / draft  # m
# gmt = kb + bmt - kg  # m
# gml = kb + bml - kg  # m

# g = 9.81 # m/s**2
# rho = 1000 # kg/m**3
# disp = wi * th * draft # m**3
# m = disp * rho # kg
# c33 = rho * g * awl # N/m
# c44 = m * g * gmt # Nm/rad
# c55 = m * g * gml # Nm/rad

# stiffnes_matrix = body.hydrostatic_stiffness.values

# print('C33: {:.2f} vs {:.2f}, ratio {:.2f}'.format(stiffnes_matrix[2,2],c33,stiffnes_matrix[2,2]/c33))
# print('C44: {:.2f} vs {:.2f}, ratio {:.2f}'.format(stiffnes_matrix[3,3],c44,stiffnes_matrix[3,3]/c44))
# print('C55: {:.2f} vs {:.2f}, ratio {:.2f}'.format(stiffnes_matrix[4,4],c55,stiffnes_matrix[4,4]/c55))
# #print(body.dofs.keys())
# body.show_matplotlib()

# solving hydrodynamics
solver = cpt.BEMSolver()
diff_prob = cpt.DiffractionProblem(body=body, wave_direction=B, water_depth=depth,omega=w)
diff_result = solver.solve(diff_prob,keep_details=(True))
rad_prob = [
    cpt.RadiationProblem(body=body, radiating_dof=dof, water_depth=depth,omega=w)
    for dof in body.dofs
    ]
rad_result = solver.solve_all(rad_prob,keep_details=(True))
#print(rad_result)

dataset = cpt.assemble_dataset(rad_result + [diff_result])
RAO = cpt.post_pro.rao(dataset, wave_direction=B, dissipation=None, stiffness=None)
pitch_RAO = np.array(np.abs(RAO.values))            # this is essentially the true pitch amplitude
# print('rao magnitude',pitch_RAO)
# print(body.dofs["Pitch"])

# # hydrodynamic and hydrostatic coefficients
# damp = dataset['radiation_damping'].sel(radiating_dof=['Pitch'],
#                                      influenced_dof=['Pitch'])
# add = dataset['added_mass'].sel(radiating_dof=['Pitch'],
#                                      influenced_dof=['Pitch'])
# stiff = dataset['hydrostatic_stiffness'].sel(radiating_dof=['Pitch'],
#                                     influenced_dof=['Pitch'])
# mass = dataset['inertia_matrix'].sel(radiating_dof=['Pitch'],
#                                      influenced_dof=['Pitch'])

# print('damping',damp)
# print('added mass', add)
# print('stiffness', stiff)
# print('mass',mass)

# Amplitude of motion of panel with largest motion
# print(np.max(np.linalg.norm(body.dofs["Surge"], axis=-1)))
# 1.0
# print(np.max(np.linalg.norm(body.dofs["Pitch"], axis=-1)))
# 11.24 --> cog, interesting

# defining the computational grid and preparing post-process data
x1 = int(-2*lambda_wave)
x2 = int(2*lambda_wave)
nx = int(lambda_wave)
y1 = int(-2*lambda_wave)
y2 = int(2*lambda_wave)
ny = int(lambda_wave)
grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
diffraction = solver.compute_free_surface_elevation(grid, diff_result)  # wave el due to diffraction
radiation = (solver.compute_free_surface_elevation(grid, rad_result[0]))*pitch_RAO  # wave el due to pitch radiation
incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)     # incident wave el
total = diffraction + radiation + incoming_fse                          # total wave el
kd = total/incoming_fse
np.savetxt('kd_050.csv', kd, delimiter=',')
np.savetxt('dif_plus_inc_data_bw.csv', total, delimiter=",")
np.savetxt('incoming_fse_data_bw.csv',incoming_fse,delimiter=',')

# plots
Z = np.real(radiation)
X = grid[0]
Y = grid[1]
plt.pcolormesh(X, Y, Z, cmap='inferno')
plt.xlabel("x")
plt.ylabel("y")
colorbar = plt.colorbar()
colorbar.set_label('Elevation')
#colorbar.set_ticks([-1, 0, 1])  # Set custom ticks
plt.tight_layout()
plt.savefig(f'elev_disturb/050elev.pdf')
plt.show()
#np.savetxt('elevation_data_rm3.csv', Z, delimiter=",")
