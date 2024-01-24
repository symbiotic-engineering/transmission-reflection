# testing hydrostatic calcs over different frequencies
import capytaine as cpt
from scipy.linalg import block_diag
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.patches as patches
from scipy.linalg import block_diag
from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity, froude_krylov_force
from capytaine.bem.airy_waves import airy_waves_free_surface_elevation

# initializing parameters
r = 1.75     # radius [m]
l = 29
x = -50
y = 0
nr = 20         # number of panels along radius
ntheta = 20     # number of panels in theta direction
nz = 20         # number of panels in z-direction
z = -r
w = 1.047        # frequency [rad/s]
B = np.pi      # wave direction [rad]
depth = 40      # average water depth at southfork

# defining mesh
body = cpt.FloatingBody(mesh=cpt.mesh_horizontal_cylinder(length=l, radius=r,center=(x,y,z),resolution=(nr,ntheta,nz),name='cyl'))
body.keep_immersed_part()
body.center_of_mass=(0,0,-1)
dofs = body.add_all_rigid_body_dofs()
body.inertia_matrix = body.compute_rigid_body_inertia()
body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()

# create array
array = body.assemble_regular_array(30,(4,1))
array.add_all_rigid_body_dofs()
array.keep_only_dofs(dofs=['0_0__Pitch','1_0__Pitch','2_0__Pitch','3_0__Pitch'])
#array.show()

# solving hydrodynamics
solver = cpt.BEMSolver()
diff_prob = cpt.DiffractionProblem(body=array, wave_direction=B, water_depth=depth,omega=w)
diff_result = solver.solve(diff_prob,keep_details=(True))
rad_prob = [
    cpt.RadiationProblem(body=array, radiating_dof=dof, water_depth=depth,omega=w)
    for dof in array.dofs
    ]
rad_result = solver.solve_all(rad_prob,keep_details=(True))
dataset = cpt.assemble_dataset(rad_result + [diff_result])

# # hydro coeffs
# damp = dataset['radiation_damping'].sel(radiating_dof=['0_0__Pitch','1_0__Pitch','2_0__Pitch','3_0__Pitch'],
#                                      influenced_dof=['0_0__Pitch','1_0__Pitch','2_0__Pitch','3_0__Pitch'])
# add = dataset['added_mass'].sel(radiating_dof=['0_0__Pitch','1_0__Pitch','2_0__Pitch','3_0__Pitch'],
#                                      influenced_dof=['0_0__Pitch','1_0__Pitch','2_0__Pitch','3_0__Pitch'])
# stiff = dataset['hydrostatic_stiffness'].sel(radiating_dof=['0_0__Pitch','1_0__Pitch','2_0__Pitch','3_0__Pitch'],
#                                      influenced_dof=['0_0__Pitch','1_0__Pitch','2_0__Pitch','3_0__Pitch'])
# mass = dataset['inertia_matrix'].sel(radiating_dof=['0_0__Pitch','1_0__Pitch','2_0__Pitch','3_0__Pitch'],
#                                      influenced_dof=['0_0__Pitch','1_0__Pitch','2_0__Pitch','3_0__Pitch'])

# post-processing
# creating mesh of free surface
x1 = -200
x2 = 200
y1 = -200
y2 = 200
nx = 100
ny = 100
grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
diffraction = solver.compute_free_surface_elevation(grid, diff_result)
radiation = sum(solver.compute_free_surface_elevation(grid, rad) for rad in rad_result[:4])
incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)
rad_dif = radiation + diffraction
total = radiation + diffraction + incoming_fse

# plots
Z = np.real(total)
X = grid[0]
Y = grid[1]
plt.pcolormesh(X, Y, Z)
plt.xlabel("x")
plt.ylabel("y")
colorbar = plt.colorbar()
colorbar.set_label('elevation')
colorbar.set_ticks([-1, 0, 1])  # Set custom ticks
plt.tight_layout()
plt.savefig(f'plots/pelamis_tot_el.pdf')
plt.show()