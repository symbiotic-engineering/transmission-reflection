# testing hydrostatic calcs over different frequencies
import capytaine as cpt
from scipy.linalg import block_diag
from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.patches as patches
from scipy.interpolate import griddata
from scipy.linalg import block_diag
from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity, froude_krylov_force
from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
import pandas as pd
import csv

# initializing parameters
wi = 40        #20         # width of box [m]
th = 10       #5          # thickness of box [m]
h = 40        #2          # height of box [m]
x = 0
y = 0
z = -18 # box center [m]
nw = 40         # number of panels along width (x)
nt = 10          # number of panels along thickness (y)
nh = 40         # number of panels along height (z)
w = 1.047        # frequency [rad/s]
B = np.pi/2     #3*np.pi/4      # wave direction [rad]
# depth = 40      # average water depth at southfork

# defining mesh
body = cpt.FloatingBody(cpt.meshes.predefined.rectangles.mesh_parallelepiped(size=(wi, th, h), resolution=(nw, nt, nh), center=(x, y, z), name='rect'))
body.keep_immersed_part()
body.center_of_mass=(0,0,-h/2)
dofs = body.add_all_rigid_body_dofs()
body.inertia_matrix = body.compute_rigid_body_inertia()
body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
#body.show_matplotlib(normal_vectors = True)

# solving hydrodynamics
solver = cpt.BEMSolver()
diff_prob = cpt.DiffractionProblem(body=body, wave_direction=B, omega=w)    # water_depth = depth
diff_result = solver.solve(diff_prob,keep_details=(True))
rad_prob = [
    cpt.RadiationProblem(body=body, radiating_dof=dof, omega=w)             # water_depth = depth
    for dof in body.dofs
    ]
rad_result = solver.solve_all(rad_prob,keep_details=(True))

dataset = cpt.assemble_dataset(rad_result + [diff_result])

# Read mesh properties and get incoming potential
# faces_centers = body.mesh.faces_centers
# phi_inc = airy_waves_potential(grid, diff_prob)
# np.savetxt('airy_phi_inc.csv', phi_inc, delimiter=",")

x1 = -100
x2 = 100
nx = 150
y1 = -100
y2 = 100
ny = 150
free = cpt.post_pro.free_surfaces.FreeSurface(x_range=(x1, x2), nx=nx, y_range=(y1, y2), ny=ny, name='test_surface')
faces_centers = free.mesh.faces_centers
phi_inc1 = airy_waves_potential(faces_centers, diff_prob)
phi_inc = phi_inc1.reshape(nx,ny)
phi_dif1 = solver.get_potential_on_mesh(diff_result, free.mesh)  #, chunk_size=50)
phi_dif = phi_dif1.reshape(nx,ny)
phi_total = phi_inc + phi_dif
grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
Z = np.real(phi_total)                #np.real((incoming_fse + fse)/incoming_fse)
X = grid[0]
Y = grid[1]
plt.pcolormesh(X, Y, Z)
plt.xlabel("x")
plt.ylabel("y")
colorbar = plt.colorbar()
colorbar.set_label('total potential')
colorbar.set_ticks([-1, 0, 1])  # Set custom ticks
plt.tight_layout()
plt.savefig('total_pot.pdf')
plt.show()

# np.savetxt('potential_data_bw.csv', phi_dif, delimiter=",")
# np.savetxt('airy_phi_inc.csv', phi_inc, delimiter=",")

# # from capytaine cookbook
# x1 = -100
# x2 = 100
# nx = 150
# y1 = -100
# y2 = 100
# ny = 150
# grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
# fse = solver.compute_free_surface_elevation(grid, diff_result)
# incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)
# total = fse + incoming_fse
# # np.savetxt('dif_plus_inc_data_bw.csv', total, delimiter=",")
# # np.savetxt('dif_el_data_bw.csv',fse,delimiter=',')
# # np.savetxt('incoming_fse_data_bw.csv',incoming_fse,delimiter=',')

# Z = np.real(total)                #np.real((incoming_fse + fse)/incoming_fse)
# X = grid[0]
# Y = grid[1]
# plt.pcolormesh(X, Y, Z)
# plt.xlabel("x")
# plt.ylabel("y")
# colorbar = plt.colorbar()
# colorbar.set_label('incoming_fse + fse')
# colorbar.set_ticks([-1, 0, 1])  # Set custom ticks
# plt.tight_layout()
# plt.savefig('breakwater_wave_transimission.pdf')
# plt.show()



# print('lambda', diff_result.wavelength)

# # damping
# damp = dataset['radiation_damping'].sel(radiating_dof=['Surge'],
#                                      influenced_dof=['Surge'])
# print('damping',damp)
# # added mass
# add = dataset['added_mass'].sel(radiating_dof=['Surge'],
#                                      influenced_dof=['Surge'])
# print('added mass',add)
# # hydrostatic stiffness
# stiff = dataset['hydrostatic_stiffness'].sel(radiating_dof=['Surge'],
#                                      influenced_dof=['Surge'])
# print('stiffness',stiff)
# # inertia matrix
# mass = dataset['inertia_matrix'].sel(radiating_dof=['Surge'],
#                                      influenced_dof=['Surge'])
# print('inertia matrix',mass)


# # Define the rectangle coordinates and size
# rectangle_x = -5
# rectangle_y = -20
# rectangle_width = 10
# rectangle_height = 40

# # Create a hatched rectangle
# rectangle = patches.Rectangle((rectangle_x, rectangle_y),
#                                rectangle_width, rectangle_height,
#                                linewidth=1, edgecolor='white', hatch='/',
#                                facecolor='none')

# # Add the rectangle to the plot
# plt.gca().add_patch(rectangle)

# cbar = fig.colorbar(CS)
# cbar.set_label('$K_{t/r}$',fontsize=17,rotation=270,labelpad=17)

# ax.set_xlabel('x [m]',fontsize=17)
# ax.set_ylabel('y [m]',fontsize=17)
# cbar.ax.tick_params(labelsize=15)
# plt.xticks(fontsize=14, rotation=90)
# plt.yticks(fontsize=14, rotation=90)