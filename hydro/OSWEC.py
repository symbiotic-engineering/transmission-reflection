def lpf(w,res):
    # hydro
    import capytaine as cpt
    from scipy.linalg import block_diag
    import numpy as np
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
    nw = 15         # number of panels along width (x)
    nt = 3          # number of panels along thickness (y)
    nh = 15         # number of panels along height (z)
    B = 0    #3*np.pi/4      # wave direction [rad]
    depth = 40      # average water depth at southfork
    g = 9.81            # gravitational constant (m/s^2)
    k = w**2/g      # wave number infinite depth (rad^2/m)
    lam = int(2*np.pi/k)    # wavelength infinite depth (m)

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


    # defining the computational grid and preparing post-process data
    x1, x2, nx, y1, y2, ny = -res*lam, res*lam, res*lam, -res*lam, res*lam, res*lam
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    diffraction = solver.compute_free_surface_elevation(grid, diff_result)  # wave el due to diffraction
    radiation = (solver.compute_free_surface_elevation(grid, rad_result[0]))*pitch_RAO  # wave el due to pitch radiation
    incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)     # incident wave el
    total = diffraction + radiation + incoming_fse                          # total wave el
    kd = total/incoming_fse
    return kd, total, incoming_fse, lam

def hydro(dataset, body):
    # hydro coeffs
    damp = (dataset['radiation_damping'].sel(radiating_dof=['Pitch'],
                                            influenced_dof=['Pitch'])).values
    add = (dataset['added_mass'].sel(radiating_dof=['Pitch'],
                                    influenced_dof=['Pitch'])).values
    add = add[0][0][0]
    damp = damp[0][0][0]
    stiff = body.hydrostatic_stiffness.values
    mass = body.inertia_matrix.values
    return damp, add, stiff, mass

