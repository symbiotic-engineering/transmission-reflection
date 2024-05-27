def lpf(w,res,xtrans,ytrans,farm):
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
                                                                                name='rect1'))
    OSWEC.keep_immersed_part()
    OSWEC.center_of_mass=(0,0,cog)
    OSWEC.add_all_rigid_body_dofs()
    OSWEC.inertia_matrix = OSWEC.compute_rigid_body_inertia()
    OSWEC.hydrostatic_stiffness = OSWEC.compute_hydrostatic_stiffness()
    OSWEC.keep_only_dofs(dofs='Pitch')

    # create array
    array = PA + OSWEC + OSWEC.translated((xtrans[0],2*ytrans[0],0),name='rect2')
    array.add_all_rigid_body_dofs()
    array.keep_only_dofs(dofs=['cyl__Heave','rect1__Pitch','rect2__Pitch'])
    # array.show_matplotlib()

    return array

# import solve
# import numpy as np
# import PTO

# xtrans = np.array([0,0])                        # x translation of bodies if farm
# ytrans = np.array([50,-50])                     # y translation of bodies if farm
# res = 2
# farm = True
# controls = True
# rad = True
# w = 1.047
# B = 0
# depth = 40

# array = lpf(w,res,xtrans,ytrans,farm)
# diff_result,rad_result,RAO_vals,lam = solve.hydro(array,B,depth,w,farm,controls)
# kd, total, incoming_fse = solve.elevation(res,lam,diff_result,rad_result,RAO_vals,farm,rad)