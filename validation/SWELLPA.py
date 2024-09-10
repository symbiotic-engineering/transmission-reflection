def lpf(w,xtrans,ytrans,depth,farm):
    # hydro
    import capytaine as cpt
    from scipy.linalg import block_diag
    import numpy as np
    from scipy.interpolate import griddata
    from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity, froude_krylov_force
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
    import pandas as pd
    import csv
    import matplotlib.pyplot as plt
    import xarray as xr

    # initializing parameters
    r = (0.256/2)     # radius [m]
    l = (0.2)
    draft = 0.110
    x = 6.79
    y = -1*4.2
    nr = 20         # number of panels along radius
    ntheta = 20     # number of panels in theta direction
    nz = 10         # number of panels in z-direction
    z = 0.5*l-draft
    B = 3*np.pi/2          # wave direction [rad]
    g = 9.81        # gravitational constant (m/s^2)
    k = w**2/g      # wave number infinite depth (rad^2/m)
    lam = int(2*np.pi/k)
    cob = -0.021
    CGz = 0.094
    CGy = 0.022

    gauge_x = np.array([6.79, 6.79, 6.79, 6.985, 6.595, 6.205, 5.815, 5.425, 6.79, 6.4, 
                        6.01, 5.62, 6.4, 10.335])
    gauge_y = -1*np.array([2.845, 3.095, 3.295, 4.2, 4.2, 4.2, 4.2, 4.2, 4.608, 4.608, 4.608, 
                        4.608, 5.87, 4.2])
    # defining mesh
    body = cpt.FloatingBody(mesh=cpt.mesh_vertical_cylinder(length=l, radius=r,center=(x,y,z),
                                                            resolution=(nr,ntheta,nz),name='cyl'),
                                                            dofs=cpt.rigid_body_dofs(rotation_center=(0,CGy,CGz)))
    body.keep_immersed_part()
    body.center_of_mass=(0,CGy,CGz)
    body.add_all_rigid_body_dofs()

    m_arm = 1.157                  # mass of arm
    m_float = 4                    # mass of float
    m = m_arm + m_float
    def mass_matrix(m,r,x=0.,y=0.):
        M = np.eye(6)*m
        M[3:,3:] *= r**2
        
        # Rotation about global origin. Adjust mass matrix accordingly:
        M[3,3] += m*x**2 
        M[4,4] += m*y**2 
        return M
    M = mass_matrix(m,r)
    body.inertia_matrix = xr.DataArray(M, dims=['influenced_dof', 'radiating_dof'],
                                        coords={'influenced_dof': list(body.dofs),'radiating_dof': list(body.dofs)},
                                        name='inertia_matrix')
    
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    body.keep_only_dofs(dofs='Heave')

    # create array
    array = body + body.translated((xtrans[0],ytrans[0],0),name='2') + body.translated((xtrans[1],ytrans[1],0),name='3')
    array.add_all_rigid_body_dofs()
    array.keep_only_dofs(dofs=['cyl__Heave','2__Heave','3__Heave'])

    if farm == False:
        array = body

    # solving hydrodynamics
    solver = cpt.BEMSolver()
    diff_prob = cpt.DiffractionProblem(body=array, water_depth=depth, wave_direction=B, omega=w)
    diff_result = solver.solve(diff_prob,keep_details=(True))
    rad_prob = [
    cpt.RadiationProblem(body=array, water_depth=depth, radiating_dof=dof, omega=w)
    for dof in array.dofs
    ]
    rad_result = solver.solve_all(rad_prob,keep_details=(True))
    dataset = cpt.assemble_dataset(rad_result + [diff_result])

    # extract hydro coeffs
    # to include off-diagonals
    A = np.squeeze(np.array([[dataset['added_mass'].sel(radiating_dof=effecting, influenced_dof=effected) for effecting in array.dofs] for effected in array.dofs]))
    B = np.squeeze(np.array([[dataset['radiation_damping'].sel(radiating_dof=effecting, influenced_dof=effected) for effecting in array.dofs] for effected in array.dofs]))
    K = array.hydrostatic_stiffness.values
    M = array.inertia_matrix.values

    # extract forces and compute exciting force
    FK = np.array([froude_krylov_force(diff_prob)[dof] for dof in array.dofs])
    dif = np.array([diff_result.forces[dof] for dof in array.dofs])
    ex_force = FK + dif
    #print('scaled ex_force',abs(0.05*ex_force))

    inertia = M + A 
    resistance = B 
    reactance = K 
    H = -(w**2)*inertia - 1j*w*resistance + reactance 

    if farm:
        RAO_controlled = np.linalg.solve(H,ex_force).ravel()
    else:
        RAO_controlled = ex_force/H
    print('RAO_controlled',RAO_controlled)

    if farm==False:
        excitationForce_SWELL = np.array([17.0876])             # pulled from their data on MATLAB
    else:                                                       # pulled from their data on MATLAB
        excitationForce_SWELL = np.array([16.7575566115100,15.2911835464861,17.2644688499461])

    # ex_force multiplied by 0.05 to match SWELL incident wave height
    FexError = (abs(0.05*ex_force) - excitationForce_SWELL)/excitationForce_SWELL
    print('F_ex error',100*FexError)

    # generating wave height and disturbance datasets (5, 11, 2.5, 6)
    x1, x2, nx, y1, y2, ny = 5, 12, 100, -1*2, -1*7, 100
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    diffraction = solver.compute_free_surface_elevation(grid, diff_result)
    multiplications = []
    if farm == False:
        for i in range(1):
            mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * RAO_controlled[i]
            multiplications.append(mult_result)
        radiation = sum(multiplications)
    else:
        for i in range(3):
            mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * RAO_controlled[i]
            multiplications.append(mult_result)
        radiation = sum(multiplications)
    incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)
    total = (incoming_fse + diffraction + radiation)*0.05       # *0.05 to match SWELL incident wave height

    import matplotlib.patheffects as path_effects
    # plots
    Z = np.real(total)
    X = grid[0]
    Y = grid[1]
    plt.pcolormesh(X, Y, Z) 
    plt.xlabel("x")
    plt.ylabel("y")
    colorbar = plt.colorbar()
    colorbar.set_label(r'Wave Elevation')
    plt.scatter(xtrans + x,ytrans + y, marker = 'o', color = 'black', s = 100)
    plt.scatter(x,y, marker = 'o', color = 'black', s = 100)
    #plt.arrow(-50, 50, 20, 0, color='black', width=0.2, head_width=5, head_length=5)
    #text = plt.text(-60, 40, 'Incident Waves', color='black', fontsize=12, ha='center', va='center')
    #text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
    plt.figure(figsize=(10, 8)) 
    plt.scatter(gauge_x,gauge_y, marker = 'o', color = 'red', s = 35)
    # x_offset = 0.05  # Adjust as needed
    y_offset = 0.05  # Adjust as needed
    #for i, (x, y) in enumerate(zip(gauge_x, gauge_y), start=1):
    #    plt.text(x, y - y_offset, f'{i}', fontsize=12, ha='center', va='top')
    # Add labels, alternating between above and below the points with specified offset
    # for i, (x, y) in enumerate(zip(gauge_x, gauge_y)):
    #     if i % 2 == 0:
    #         plt.text(x, y - y_offset, f'({x:.2f}, {y:.2f})', fontsize=10, ha='center', va='top')
    #     else:
    #         plt.text(x, y + y_offset, f'({x:.2f}, {y:.2f})', fontsize=10, ha='center', va='bottom')
    plt.xlabel('X [m]',fontsize='14')
    plt.ylabel('Y [m]',fontsize='14')
    plt.xticks(fontsize='12')
    plt.yticks(fontsize='12')
    plt.tight_layout()
    plt.savefig('validation_sensors.pdf')
    #plt.show()

    points = np.column_stack((gauge_x, gauge_y))
    # Interpolate 'total' onto the gauge points
    elevation_at_gauges = griddata((X.ravel(), Y.ravel()), total.ravel(), points, method='linear')
    print('wave elevation',np.abs(elevation_at_gauges))
    print('wave amplitude',np.abs(elevation_at_gauges)/2)
    # total_at_gauges now contains the values of 'total' at the gauge locations specified by gauge_x and gauge_y

    return total, incoming_fse, lam, elevation_at_gauges