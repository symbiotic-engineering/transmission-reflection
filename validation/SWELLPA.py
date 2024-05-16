def lpf(w,res,xtrans,ytrans,depth,froude_number,farm):
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
    r = (0.256/2)/froude_number     # radius [m]
    l = (0.2)/froude_number
    draft = 0.110/froude_number
    x = 6.79/froude_number
    y = -1*4.2/froude_number
    nr = 20         # number of panels along radius
    ntheta = 20     # number of panels in theta direction
    nz = 10         # number of panels in z-direction
    z = 0.5*l-draft
    B = 3*np.pi/2          # wave direction [rad]
    g = 9.81        # gravitational constant (m/s^2)
    k = w**2/g      # wave number infinite depth (rad^2/m)
    lam = int(2*np.pi/k)
    cob = -0.021/froude_number
    CGz = 0.094/froude_number
    CGy = 0.022/froude_number

    gauge_x = np.array([6.79, 6.79, 6.79, 6.985, 6.595, 6.205, 5.815, 5.425, 6.79, 6.4, 
                        6.01, 5.62, 6.4, 10.335])/froude_number #, 6.79, 6.4, 6.01, 5.62, 6.4])/froude_number
    gauge_y = -1*np.array([2.845, 3.095, 3.295, 4.2, 4.2, 4.2, 4.2, 4.2, 4.608, 4.608, 4.608, 
                        4.608, 5.87, 4.2])/froude_number #, 4.2, 4.2, 4.2, 4.2, 5.07])/froude_number
    # defining mesh
    body = cpt.FloatingBody(mesh=cpt.mesh_vertical_cylinder(length=l, radius=r,center=(x,y,z),
                                                            resolution=(nr,ntheta,nz),name='cyl'),
                                                            dofs=cpt.rigid_body_dofs(rotation_center=(0,CGy,CGz)))
    body.keep_immersed_part()
    body.center_of_mass=(0,CGy,CGz)
    body.add_all_rigid_body_dofs()

    m_arm = 1.157/froude_number**3                   # mass of arm
    m_float = 4/froude_number**3                     # mass of float
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
    
    # # initializing parameters
    # wi = 19.3         # width of box [m]
    # th = 1          # thickness of box [m]
    # h = 2           # height of box [m]
    # xwall = 0
    # ywall = -10
    # z = -0.5        # box center [m]
    # nw = 15         # number of panels along width (x)
    # nt = 2         # number of panels along thickness (y)
    # nh = 3          # number of panels along height (z)

    # # defining mesh
    # wall = cpt.FloatingBody(cpt.meshes.predefined.rectangles.mesh_parallelepiped(size=(th, wi, h), 
    #                                                                              resolution=(nt, nw, nh), 
    #                                                                              center=(xwall, ywall, z), 
    #                                                                              name='wall'))
    # wall.keep_immersed_part()
    # wall.center_of_mass=(0,0,-h/2)
    # wall.add_all_rigid_body_dofs()
    # wall.inertia_matrix = wall.compute_rigid_body_inertia()
    # wall.hydrostatic_stiffness = wall.compute_hydrostatic_stiffness()
    # wall_dist = 14.6 + (th/2)

    # back_wall = cpt.FloatingBody(cpt.meshes.predefined.rectangles.mesh_parallelepiped(size=(wi, th, h), 
    #                                                                              resolution=(nw, nt, nh), 
    #                                                                              center=(8, -20, z), 
    #                                                                              name='backwall'))
    # back_wall.keep_immersed_part()
    # back_wall.center_of_mass=(0,0,-h/2)
    # back_wall.add_all_rigid_body_dofs()
    # back_wall.inertia_matrix = back_wall.compute_rigid_body_inertia()
    # back_wall.hydrostatic_stiffness = back_wall.compute_hydrostatic_stiffness()

    # walls = wall + wall.translated((wall_dist,0,0),name='wall2') + back_wall
    # walls.add_all_rigid_body_dofs()
    # walls.show_matplotlib()


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
    RAO = cpt.post_pro.rao(dataset, wave_direction=B, dissipation=None, stiffness=None)
    RAO = np.array(np.abs(RAO.values))            # this is essentially the true Heave amplitude
    print('rao',RAO)

    #wall_diff = cpt.DiffractionProblem(body=walls, water_depth=depth, wave_direction=B, omega=w)
    #wall_diff_result = solver.solve(wall_diff,keep_details=(True))

    # hydro coeffs
    damp = np.array([dataset['radiation_damping'].sel(radiating_dof=dof,
                                            influenced_dof=dof) for dof in array.dofs])
    add = np.array([dataset['added_mass'].sel(radiating_dof=dof,
                                            influenced_dof=dof) for dof in array.dofs])
    stiff = body.hydrostatic_stiffness.values
    mass = body.inertia_matrix.values
    inertial_term = [-w**2*(int(mass) + A) for A in add]
    mech_term = [-1j*w*B + int(stiff) for B in damp]
    transfer_matrix = [inertial + mech for inertial, mech in zip(inertial_term, mech_term)]
    #print('H_ij', transfer_matrix)
    if farm==False:
        excitationForce_SWELL = np.array([17.3403])/(froude_number**3)
        FK = froude_krylov_force(diff_prob)['Heave']
        diff_force = diff_result.forces['Heave']
    else:
        excitationForce_SWELL = np.array([17.0947148425932, 15.7179554013022, 17.5956099088629])/(froude_number**3)
        FK1 = froude_krylov_force(diff_prob)['cyl__Heave']
        FK2 = froude_krylov_force(diff_prob)['2__Heave']
        FK3 = froude_krylov_force(diff_prob)['3__Heave']
        FK = np.array([FK1,FK2,FK3])
        diff_force = np.array([diff_result.forces['cyl__Heave'],diff_result.forces['2__Heave'],
                                diff_result.forces['3__Heave']])
    
    F_ex = (FK + diff_force)*0.05
    print(abs(F_ex))
    FexError = (abs(F_ex) - excitationForce_SWELL)/excitationForce_SWELL
    print('error',FexError)

    experimental_RAO = np.abs([force / term for force, term in zip(excitationForce_SWELL*froude_number**3, transfer_matrix)])
    print('RAO',experimental_RAO)

    # generating wave height and disturbance datasets (5, 11, 2.5, 6)
    x1, x2, nx, y1, y2, ny = 5/froude_number, 10/froude_number, 300, -1*2/froude_number, -1*15/froude_number, 300 #-res*lam, res*lam, res*lam, -res*lam, res*lam, res*lam
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    # wall_diffraction = solver.compute_free_surface_elevation(grid, wall_diff_result)
    diffraction = solver.compute_free_surface_elevation(grid, diff_result)
    multiplications = []
    if farm == False:
        for i in range(1):
            mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * experimental_RAO[i,0] #Heave_RAO[0,i]
            multiplications.append(mult_result)
        radiation = sum(multiplications)
    else:
        for i in range(3):
            mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * experimental_RAO[i,0] #Heave_RAO[0,i]
            multiplications.append(mult_result)
        radiation = sum(multiplications)
    incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)
    total = (diffraction + radiation + incoming_fse)*0.05/(froude_number)
    kd = total/incoming_fse

    # plots
    Z = np.real(total)
    X = grid[0]
    Y = grid[1]
    plt.pcolormesh(X, Y, Z) #, cmap=cmap,vmin=0,vmax=2.5)
    plt.xlabel("x")
    plt.ylabel("y")
    colorbar = plt.colorbar()
    colorbar.set_label(r'Wave Elevation')
    plt.scatter(xtrans + x,ytrans + y, marker = 'o', color = 'black', s = 100)
    plt.scatter(x,y, marker = 'o', color = 'black', s = 100)
    plt.scatter(gauge_x,gauge_y, marker = 'o', color = 'red', s = 25)
    # plt.arrow(-50, 50, 20, 0, color='black', width=0.2, head_width=5, head_length=5)
    # plt.text(-60, 40, 'Incident Waves', color='black', fontsize=12, ha='center', va='center')
    plt.tight_layout()
    plt.savefig('validation.pdf')
    plt.show()

    # After calculating 'total' in your function, add the following lines:
    # Create a grid of points from gauge_x and gauge_y
    points = np.column_stack((gauge_x, gauge_y))
    # Interpolate 'total' onto the gauge points
    elevation_at_gauges = griddata((X.ravel(), Y.ravel()), total.ravel(), points, method='linear')
    print(np.abs(elevation_at_gauges))
    # total_at_gauges now contains the values of 'total' at the gauge locations specified by gauge_x and gauge_y

    return kd, total, incoming_fse, lam, elevation_at_gauges