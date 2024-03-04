def lpf(w,res,xtrans,ytrans):
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
    nw = 10         # number of panels along width (x)
    nt = 3          # number of panels along thickness (y)
    nh = 10         # number of panels along height (z)
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
    body.add_all_rigid_body_dofs()
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    body.keep_only_dofs(dofs='Pitch')

    # create array
    array = body + body.translated((xtrans[0],ytrans[0],0),name='2') + body.translated((xtrans[1],ytrans[1],0),name='3')
    array.add_all_rigid_body_dofs()
    array.keep_only_dofs(dofs=['rect__Pitch','2__Pitch','3__Pitch'])
    # array.show_matplotlib()

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
    RAO = cpt.post_pro.rao(dataset, wave_direction=B, dissipation=None, stiffness=None)
    pitch_RAO = np.array(np.abs(RAO.values))            # this is essentially the true pitch amplitude

    # defining the computational grid and preparing post-process data
    x1, x2, nx, y1, y2, ny = -res*lam, res*lam, res*lam, -res*lam, res*lam, res*lam
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    diffraction = solver.compute_free_surface_elevation(grid, diff_result)  # wave el due to diffraction
    multiplications = []
    for i in range(3):
        mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * pitch_RAO[0,i]
        multiplications.append(mult_result)
    radiation = sum(multiplications)
    incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)     # incident wave el
    total = diffraction + radiation + incoming_fse                          # total wave el
    kd = total/incoming_fse
    return kd, total, incoming_fse, lam