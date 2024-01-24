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
nw = 40         # number of panels along width (x)
nt = 2          # number of panels along thickness (y)
nh = 40         # number of panels along height (z)
w = 1.047        # frequency [rad/s]
B = np.pi/2     #3*np.pi/4      # wave direction [rad]
depth = 40      # average water depth at southfork

# defining mesh
body = cpt.FloatingBody(cpt.meshes.predefined.rectangles.mesh_parallelepiped(size=(wi, th, h), resolution=(nw, nt, nh), center=(x, y, z), name='rect'))
body.keep_immersed_part()
body.center_of_mass=(0,0,cog)
dofs = body.add_all_rigid_body_dofs()
body.inertia_matrix = body.compute_rigid_body_inertia()
body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()

#body.show_matplotlib(normal_vectors = True)

# solving hydrodynamics
solver = cpt.BEMSolver()
diff_prob = cpt.DiffractionProblem(body=body, wave_direction=B, water_depth=depth,omega=w)
diff_result = solver.solve(diff_prob,keep_details=(True))
rad_prob = [
    cpt.RadiationProblem(body=body, radiating_dof=dof, water_depth=depth,omega=1.047)
    for dof in body.dofs
    ]
rad_result = solver.solve_all(rad_prob,keep_details=(True))

dataset = cpt.assemble_dataset(rad_result + [diff_result])

# hydrodynamic and hydrostatic coefficients
# damp = dataset['radiation_damping'].sel(radiating_dof=['Surge'],
#                                      influenced_dof=['Surge'])
# add = dataset['added_mass'].sel(radiating_dof=['Surge'],
#                                      influenced_dof=['Surge'])
# stiff = dataset['hydrostatic_stiffness'].sel(radiating_dof=['Surge'],
#                                      influenced_dof=['Surge'])
# mass = dataset['inertia_matrix'].sel(radiating_dof=['Surge'],
#                                      influenced_dof=['Surge'])

# defining the computational grid and preparing post-process data
x1, x2, nx = -100, 100, 150
y1, y2, ny = -100, 100, 150
grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
diffraction = solver.compute_free_surface_elevation(grid, diff_result)  # wave el due to diffraction
radiation = solver.compute_free_surface_elevation(grid, rad_result[0])  # wave el due to surge radiation
incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)     # incident wave el
total = diffraction + radiation + incoming_fse                          # total wave el
upeffect = diffraction + radiation                                      # reflected wave elevation

# plots
Z = np.real(total)
X = grid[0]
Y = grid[1]
plt.pcolormesh(X, Y, Z)
plt.xlabel("x")
plt.ylabel("y")
colorbar = plt.colorbar()
colorbar.set_label('fse')
colorbar.set_ticks([-1, 0, 1])  # Set custom ticks
plt.tight_layout()
plt.savefig(f'plots/tot_el.pdf')
plt.show()
#np.savetxt('elevation_data_rm3.csv', Z, delimiter=",")
