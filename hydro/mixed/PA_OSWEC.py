def lpf(w,res,xtrans,ytrans):
    # hydro
    import capytaine as cpt
    from scipy.linalg import block_diag
    import numpy as np
    from scipy.interpolate import griddata
    from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity, froude_krylov_force
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
    import pandas as pd
    import matplotlib.pyplot as plt
    import csv

    # initializing parameters
    # point absorber
    r = 10     # radius [m]
    l = 5
    x = 0
    y = 0
    nr = 20         # number of panels along radius
    ntheta = 10     # number of panels in theta direction
    nz = 5         # number of panels in z-direction
    z = 0

    # oscillating surge
    wi = 25         # width of flap [m]
    th = 1          # thickness of flap [m]
    h = 19          # height of flap [m], draft = 16 m
    draft = 16 # m
    z_flap = 0.5*h-draft # box center [m]
    cog = -11.4     # center of gravity [m]
    nw = 10         # number of panels along width (x)
    nt = 3          # number of panels along thickness (y)
    nh = 10         # number of panels along height (z)

    # wave parameters
    B = 0           # wave direction [rad]
    g = 9.81        # gravitational constant (m/s^2)
    k = w**2/g      # wave number infinite depth (rad^2/m)
    lam = int(2*np.pi/k)

    # defining point absorber mesh
    PA = cpt.FloatingBody(mesh=cpt.mesh_vertical_cylinder(length=l, radius=r,center=(x,y,z),
                                                            resolution=(nr,ntheta,nz),name='cyl'))
    PA.keep_immersed_part()
    PA.center_of_mass=(0,0,z)
    PA.add_all_rigid_body_dofs()
    PA.inertia_matrix = PA.compute_rigid_body_inertia()
    PA.hydrostatic_stiffness = PA.compute_hydrostatic_stiffness()
    PA.keep_only_dofs(dofs='Heave')

    # defining oscillating surge mesh
    OSWEC = cpt.FloatingBody(cpt.meshes.predefined.rectangles.mesh_parallelepiped(size=(th, wi, h), 
                                                                                resolution=(nt, nw, nh), 
                                                                                center=(xtrans[1], ytrans[1], z),
                                                                                name='rect'))
    OSWEC.keep_immersed_part()
    OSWEC.center_of_mass=(0,0,cog)
    OSWEC.add_all_rigid_body_dofs()
    OSWEC.inertia_matrix = OSWEC.compute_rigid_body_inertia()
    OSWEC.hydrostatic_stiffness = OSWEC.compute_hydrostatic_stiffness()
    OSWEC.keep_only_dofs(dofs='Pitch')

    # create array
    array = PA + OSWEC
    array.add_all_rigid_body_dofs()
    array.keep_only_dofs(dofs=['cyl__Heave','rect__Pitch'])
    #array.show_matplotlib()

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
    body_RAO = np.array(np.abs(RAO.values))            # this is essentially the true pitch amplitude
    print('heave_RAO',body_RAO[0,0])

    # generating wave height and disturbance datasets
    x1, x2, nx, y1, y2, ny = -res*lam, res*lam, res*lam, -res*lam, res*lam, res*lam
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    diffraction = solver.compute_free_surface_elevation(grid, diff_result)
    radiation_PA = solver.compute_free_surface_elevation(grid, rad_result[0]) * body_RAO[0,0]
    radiation_OSWEC = solver.compute_free_surface_elevation(grid, rad_result[1]) * body_RAO[0,1]
    radiation = radiation_PA + radiation_OSWEC
    incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)
    total = diffraction + radiation + incoming_fse
    kd = total/incoming_fse

    Z = np.real(radiation)
    X = grid[0]
    Y = grid[1]
    plt.pcolormesh(X, Y, Z)
    plt.xlabel("x")
    plt.ylabel("y")
    colorbar = plt.colorbar()
    colorbar.set_label('Radiated Wave Field')
    plt.tight_layout()
    plt.savefig('radiation.pdf')
    plt.show()
    return kd, total, incoming_fse, lam