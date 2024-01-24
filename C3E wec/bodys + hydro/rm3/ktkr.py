import numpy as np
import capytaine as cpt
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.linalg import block_diag
from scipy.interpolate import griddata
from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity, froude_krylov_force
from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
import pandas as pd
import csv

def wave_num(w):
    # used to find infinite depth wave number (k) and wavelength (lam)
    g = 9.81
    k = w**2 / g
    lam = 2 * np.pi / k
    return [k, lam]

def hydrodynamic_calculations(w, B, depth, r, l, x, y, nr, ntheta, nz, z):
    # used to create body mesh, solve the radiation and diffraction problems,
    # and extract wave elevation in front and in lee of the body

    # defining the mesh
    body = cpt.FloatingBody(mesh=cpt.mesh_vertical_cylinder(length=l, radius=r, center=(x, y, z),
                                                            resolution=(nr, ntheta, nz), name='cyl'))
    body.keep_immersed_part()
    body.center_of_mass = (0, 0, -l/2)
    body.add_all_rigid_body_dofs()
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    
    # solving diffraction and radiation problems
    solver = cpt.BEMSolver()
    diff_prob = cpt.DiffractionProblem(body=body, wave_direction=B, water_depth=depth, omega=w)
    diff_result = solver.solve(diff_prob, keep_details=(True))
    rad_prob = [cpt.RadiationProblem(body=body, radiating_dof=dof, water_depth=depth, omega=w)
                for dof in body.dofs]
    rad_result = solver.solve_all(rad_prob, keep_details=(True))
    dataset = cpt.assemble_dataset(rad_result + [diff_result])

    # defining the computational grid and preparing post-process data
    x1, x2, nx = -100, 100, 150
    y1, y2, ny = -100, 100, 150
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    diffraction = solver.compute_free_surface_elevation(grid, diff_result)  # wave el due to diffraction
    radiation = solver.compute_free_surface_elevation(grid, rad_result[2])  # wave el due to heave radiation
    incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)     # incident wave el
    total = diffraction + radiation + incoming_fse                          # total wave el
    upeffect = diffraction + radiation                                      # reflected wave elevation

    return incoming_fse, total, upeffect

def post_process_wave_data(xR, xD, upeffect, incoming_fse, total):
    # extracting the center line of wave elevation for computation
    z_up = np.real(upeffect[74,:75])            # just the rad + dif upstream
    z_down = np.real(total[74,75:])             # total wave elevation
    zinc_up = np.real(incoming_fse[74,:75])     # incident wave elevation
    zinc_down = np.real(incoming_fse[74,75:])   # incident wave elevation

    # Dictionary to store variables and their wave_height
    variables = {'z_up': z_up, 'z_down': z_down, 'zinc_up': zinc_up, 'zinc_down': zinc_down}
    wave_height = {}
    for variable_name, variable in variables.items():
        # Your logic to determine x_array and calculate wave_height...
        x_array = xR if variable_name.endswith('_up') else xD if variable_name.endswith('_down') else None
        if x_array is not None:
            indices = np.where((x_array < -25) | (x_array > 25))[0]
            if len(indices) > 0:
                peak = np.max(variable[indices])
                trough = np.min(variable[indices])
                difference = peak - trough
                # Store the difference in the dictionary
                wave_height[variable_name] = difference

    # find the transmission and reflection coefficients
    Kt = Kr = None  # Assign default values
    if 'z_up' in wave_height and 'zinc_up' in wave_height and wave_height['zinc_up'] != 0:
        Kr = wave_height['z_up'] / wave_height['zinc_up']
    else:
        print('Division cannot be performed.')
    if 'z_down' in wave_height and 'zinc_down' in wave_height and wave_height['zinc_down'] != 0:
        Kt = wave_height['z_down'] / wave_height['zinc_down']
    else:
        print('Division for z_down and zinc_down cannot be performed.')

    return Kt, Kr

# Define parameters
w = np.array([0.5, 0.75, 0.85, 0.9, 1.047, 1.2, 1.35, 1.5, 1.75])   # wave frequency
r = 10              # radius of body
l = 5               # length of body
x = 0               # x-position of body
y = 0               # y-position of body 
z = 0               # z-position of body
nr = 40             # number of panels in r-dimension
ntheta = 15         # number of panels in theta-dimension
nz = 10             # number of panels in z-dimension
B = 0               # wave heading --> note to olivia: investigate this
depth = 40          # ocean depth (not used in this code yet, only for finite depth)
xR = np.linspace(-100, 0, 75)       # upstream x-values
xD = np.linspace(0, 100, 75)        # downstream x-values
# Initialize empty lists for Kt and Kr
Kt_values = []
Kr_values = []

# Loop through frequencies
for freq in w:
    # Calculate wave numbers and lengths
    [k, lam] = wave_num(freq)
    # Perform hydrodynamic calculations
    incoming_fse, total, upeffect = hydrodynamic_calculations(freq, B, depth, r, l, x, y, nr, ntheta, nz, z)
    # Calculate Kt and Kr for the current frequency
    Kt, Kr = post_process_wave_data(xR, xD, upeffect, incoming_fse, total)
    # Append the values to the lists
    Kt_values.append(Kt)
    Kr_values.append(Kr)

Kt = np.array(Kt_values)
Kr = np.array(Kr_values)

plt.plot(w, Kt, color='blue', label='Transmission Coefficient')
plt.plot(w, Kr, color='red', label='Reflection Coefficient')
plt.legend()
plt.xlabel('w [rad/s]')
plt.ylabel('Coefficient')
plt.savefig('tran_ref_coeff.pdf')
plt.show()