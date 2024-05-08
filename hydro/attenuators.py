def lpf(w,res,N,D,farm):
    # hydro
    # N = number of bodies = 4
    # D = distance btwn bodies = 30
    import capytaine as cpt
    from scipy.linalg import block_diag
    import numpy as np
    from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity, froude_krylov_force
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation

    # initializing parameters
    r = 1.75     # radius [m]
    l = 29
    x = -50
    y = 0
    nr = 2        # number of panels along radius
    ntheta = 5     # number of panels in theta direction
    nz = 10         # number of panels in z-direction
    z = 0
    B = 0      # wave direction [rad]
    depth = 40      # average water depth at southfork
    g = 9.81            # gravitational constant (m/s^2)
    k = w**2/g      # wave number infinite depth (rad^2/m)
    lam = int(2*np.pi/k)    # wavelength infinite depth (m)

    # defining mesh
    body = cpt.FloatingBody(mesh=cpt.mesh_horizontal_cylinder(length=l, radius=r,center=(x,y,z),
                                                            resolution=(nr,ntheta,nz),name='cyl'))
    body.keep_immersed_part()
    body.center_of_mass=(0,0,0)
    body.add_all_rigid_body_dofs()
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()

    # create whole body
    array = body.assemble_regular_array(D,(4,N))
    array.add_all_rigid_body_dofs()
    array.keep_only_dofs(dofs=['0_0__Pitch','1_0__Pitch','2_0__Pitch','3_0__Pitch','0_1__Pitch','1_1__Pitch',
                               '2_1__Pitch','3_1__Pitch','0_2__Pitch','1_2__Pitch','2_2__Pitch','3_2__Pitch'])

    if farm == False:
        array = body.assemble_regular_array(30,(4,1))
        array.add_all_rigid_body_dofs()
        array.keep_only_dofs(dofs=['0_0__Pitch','1_0__Pitch','2_0__Pitch','3_0__Pitch'])
    
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
    pitch_RAO = np.abs(RAO.values)

    # solving for wave height and disturbance datasets
    x1, x2, nx, y1, y2, ny = -res*lam, res*lam, res*lam, -res*lam, res*lam, res*lam
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    diffraction = solver.compute_free_surface_elevation(grid, diff_result)
    multiplications = []
    if farm == False:
        for i in range(4):
            mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * pitch_RAO[0,i]
            multiplications.append(mult_result)
    else:
        for i in range(4*N):
            mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * pitch_RAO[0,i]
            multiplications.append(mult_result)
    radiation = sum(multiplications)
    incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)
    total = radiation + diffraction + incoming_fse
    kd = total/incoming_fse

    return kd, total, incoming_fse, lam