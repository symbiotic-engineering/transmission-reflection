# testing hydrostatic calcs over different frequencies
import capytaine as cpt
from scipy.linalg import block_diag
from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.patches as patches
from matplotlib.colors import Normalize
from scipy.interpolate import griddata
from scipy.linalg import block_diag
from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity, froude_krylov_force
from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
import pandas as pd
import csv

# initializing parameters
# THIS IS FOR THE LIMITING CASE i.e. pretend in flume taking up entire domain
wi = 40        #20         # width of box [m]
th = 10       #5          # thickness of box [m]
h = 40        #2          # height of box [m]
x = 0
y = 0
z = -18 # box center [m]
nw = 15         # number of panels along width (x)
nt = 5          # number of panels along thickness (y)
nh = 15         # number of panels along height (z)
w = 0.50     # frequency [rad/s]
B = np.pi/2     #3*np.pi/4      # wave direction [rad]
def wave_num(w):
    # for infinite depth
    g = 9.81            # gravitational constant (m/s^2)
    k_inf = w**2/g      # wave number infinite depth (rad^2/m)
    lam_inf = 2*np.pi/k_inf    # wavelength infinite depth (m)
    return k_inf, lam_inf

k, lambda_wave = wave_num(w)
# depth = 40      # average water depth at southfork

# defining mesh
body = cpt.FloatingBody(cpt.meshes.predefined.rectangles.mesh_parallelepiped(size=(wi, th, h), resolution=(nw, nt, nh), center=(x, y, z), name='rect'))
body.keep_immersed_part()
body.center_of_mass=(0,0,-h/2)
body.add_all_rigid_body_dofs()
body.inertia_matrix = body.compute_rigid_body_inertia()
body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
#body.show_matplotlib(normal_vectors = True)

walls = cpt.FloatingBody(cpt.meshes.predefined.rectangles.mesh_parallelepiped(size=(th, 4*wi, h), resolution=(nt, 2*nw, nh), center=((th-3)+(wi/2), y, z), name='wall'))
walls.keep_immersed_part()
walls.center_of_mass=(0,0,-h/2)
walls.add_all_rigid_body_dofs()
walls.inertia_matrix = walls.compute_rigid_body_inertia()
walls.hydrostatic_stiffness = walls.compute_hydrostatic_stiffness()

all = body + walls + walls.translated((-(2*((th-3) + (wi/2))),0,0),name='wall2')
all.show_matplotlib(normal_vectors = True)


# solving hydrodynamics
solver = cpt.BEMSolver()
diff_prob = cpt.DiffractionProblem(body=all, wave_direction=B, omega=w)    # water_depth = depth
diff_result = solver.solve(diff_prob,keep_details=(True))
rad_prob = [
    cpt.RadiationProblem(body=all, radiating_dof=dof, omega=w)             # water_depth = depth
    for dof in all.dofs
    ]
rad_result = solver.solve_all(rad_prob,keep_details=(True))

dataset = cpt.assemble_dataset(rad_result + [diff_result])

# Read mesh properties and get incoming potential
# faces_centers = body.mesh.faces_centers
# phi_inc = airy_waves_potential(grid, diff_prob)
# np.savetxt('airy_phi_inc.csv', phi_inc, delimiter=",")

# x1 = -200
# x2 = 200
# nx = 150
# y1 = -200
# y2 = 200
# ny = 150
# grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
# free = cpt.post_pro.free_surfaces.FreeSurface(x_range=(x1, x2), nx=nx, y_range=(y1, y2), ny=ny, name='test_surface')
# faces_centers = free.mesh.faces_centers
# phi_inc1 = airy_waves_potential(faces_centers, diff_prob)
# phi_inc = phi_inc1.reshape(nx,ny)
# phi_dif1 = solver.get_potential_on_mesh(diff_result,free.mesh)     #, chunk_size=50)
# phi_dif = phi_dif1.reshape(nx,ny)
# phi_total = phi_inc + phi_dif
# Z = np.real(phi_total)                #np.real((incoming_fse + fse)/incoming_fse)
# X = grid[0]
# Y = grid[1]
# plt.pcolormesh(X, Y, Z)
# plt.xlabel("x")
# plt.ylabel("y")
# colorbar = plt.colorbar()
# colorbar.set_label('total potential')
# colorbar.set_ticks([-1, 0, 1])  # Set custom ticks
# plt.tight_layout()
# plt.savefig('total_pot.pdf')
# plt.show()

# np.savetxt('potential_data_bw.csv', phi_dif, delimiter=",")
# np.savetxt('airy_phi_inc.csv', phi_inc, delimiter=",")

# from capytaine cookbook
x1 = int(-2*lambda_wave)
x2 = int(2*lambda_wave)
nx = int(lambda_wave)
y1 = int(-2*lambda_wave)
y2 = int(2*lambda_wave)
ny = int(lambda_wave)
grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
fse = solver.compute_free_surface_elevation(grid, diff_result)
incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)
total = fse + incoming_fse
kd = total/incoming_fse
np.savetxt('kd_050test.csv', kd, delimiter=',')
np.savetxt('dif_plus_inc_data_bw.csv', total, delimiter=",")
np.savetxt('dif_el_data_bw.csv',fse,delimiter=',')
np.savetxt('incoming_fse_data_bw.csv',incoming_fse,delimiter=',')

Z = np.real(total)                #np.real((incoming_fse + fse)/incoming_fse)
X = grid[0]
Y = grid[1]
#norm = Normalize(vmin=0, vmax=2)
plt.pcolormesh(X, Y, Z, cmap='inferno') #,norm=norm)
plt.xlabel("x")
plt.ylabel("y")
colorbar = plt.colorbar()
colorbar.set_label('Elevation')
#colorbar.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2])  # Set custom ticks
plt.tight_layout()
plt.savefig(f'field_plots/050test_elevation.pdf')
plt.show()