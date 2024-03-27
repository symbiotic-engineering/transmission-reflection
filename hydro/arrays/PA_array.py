def lpf(w,res,xtrans,ytrans):
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
    nr = 20         # number of panels along radius
    ntheta = 10     # number of panels in theta direction
    nz = 5         # number of panels in z-direction
    z = 0
    B = 0           # wave direction [rad]
    g = 9.81        # gravitational constant (m/s^2)
    k = w**2/g      # wave number infinite depth (rad^2/m)
    lam = int(2*np.pi/k)

    # defining mesh
    body = cpt.FloatingBody(mesh=cpt.mesh_vertical_cylinder(length=l, radius=r,center=(x,y,z),
                                                            resolution=(nr,ntheta,nz),name='cyl'))
    body.keep_immersed_part()
    body.center_of_mass=(0,0,z)
    body.add_all_rigid_body_dofs()
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    body.keep_only_dofs(dofs='Heave')

    # create array
    array = body + body.translated((xtrans[0],ytrans[0],0),name='2') + body.translated((xtrans[1],ytrans[1],0),name='3')
    array.add_all_rigid_body_dofs()
    array.keep_only_dofs(dofs=['cyl__Heave','2__Heave','3__Heave'])
    # array.show_matplotlib()

    # solving hydrodynamics
    solver = cpt.BEMSolver()
    diff_prob = cpt.DiffractionProblem(body=array, wave_direction=B, omega=w)
    diff_result = solver.solve(diff_prob,keep_details=(True))
    rad_prob = [
    cpt.RadiationProblem(body=array, radiating_dof=dof, omega=w)
    for dof in array.dofs
    ]
    rad_result = solver.solve_all(rad_prob,keep_details=(True))
    dataset = cpt.assemble_dataset(rad_result + [diff_result])
    RAO = cpt.post_pro.rao(dataset, wave_direction=B, dissipation=None, stiffness=None)
    heave_RAO = np.array(np.abs(RAO.values))            # this is essentially the true pitch amplitude

    # generating wave height and disturbance datasets
    x1, x2, nx, y1, y2, ny = -res*lam, res*lam, res*lam, -res*lam, res*lam, res*lam
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    diffraction = solver.compute_free_surface_elevation(grid, diff_result)
    multiplications = []
    for i in range(3):
        mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * heave_RAO[0,i]
        multiplications.append(mult_result)
    radiation = sum(multiplications)
    incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)
    total = diffraction + radiation + incoming_fse
    kd = total/incoming_fse
    return kd, total, incoming_fse, lam