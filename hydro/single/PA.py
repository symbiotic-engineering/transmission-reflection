def lpf(w,res):
    # hydro
    import capytaine as cpt
    from scipy.linalg import block_diag
    import numpy as np
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
    B = 0           # wave direction [rad]
    g = 9.81        # gravitational constant (m/s^2)
    k = w**2/g      # wave number infinite depth (rad^2/m)
    lam = int(2*np.pi/k)

    # defining mesh
    body = cpt.FloatingBody(mesh=cpt.mesh_vertical_cylinder(length=l, radius=r,center=(x,y,z),resolution=(nr,ntheta,nz),name='cyl'))
    body.keep_immersed_part()
    body.center_of_mass=(0,0,z)
    body.add_all_rigid_body_dofs()
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    body.keep_only_dofs(dofs='Heave')

    # solving hydrodynamics
    solver = cpt.BEMSolver()
    diff_prob = cpt.DiffractionProblem(body=body, wave_direction=B, omega=w)
    diff_result = solver.solve(diff_prob,keep_details=(True))
    rad_prob = [
    cpt.RadiationProblem(body=body, radiating_dof=dof, omega=w)
    for dof in body.dofs
    ]
    rad_result = solver.solve_all(rad_prob,keep_details=(True))
    dataset = cpt.assemble_dataset(rad_result + [diff_result])

    # generating wave height and disturbance datasets
    x1, x2, nx, y1, y2, ny = -res*lam, res*lam, res*lam, -res*lam, res*lam, res*lam
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    diffraction = solver.compute_free_surface_elevation(grid, diff_result)
    chosen_rad_result = rad_result[0]
    radiation = solver.compute_free_surface_elevation(grid, chosen_rad_result)
    incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)
    total = diffraction + radiation + incoming_fse
    kd = total/incoming_fse
    return kd, total, incoming_fse, lam

def hydro(dataset, body):
    # hydro coeffs
    damp = (dataset['radiation_damping'].sel(radiating_dof=['Heave'],
                                            influenced_dof=['Heave'])).values
    add = (dataset['added_mass'].sel(radiating_dof=['Heave'],
                                    influenced_dof=['Heave'])).values
    add = add[0][0][0]
    damp = damp[0][0][0]
    stiff = body.hydrostatic_stiffness.values
    mass = body.inertia_matrix.values
    return damp, add, stiff, mass
