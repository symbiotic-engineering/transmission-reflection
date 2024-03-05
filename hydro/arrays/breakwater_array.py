# hydrostatic calcs
def lpf(w, res,xtrans,ytrans):
    import capytaine as cpt
    from scipy.linalg import block_diag
    from scipy.spatial import Delaunay
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity, froude_krylov_force
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
    import pandas as pd
    import csv

    # initializing parameters
    wi = 20         # width of box [m]
    th = 5          # thickness of box [m]
    h = 2           # height of box [m]
    x = 0
    y = 0
    z = -0.5        # box center [m]
    nw = 10         # number of panels along width (x)
    nt = 5         # number of panels along thickness (y)
    nh = 2          # number of panels along height (z)
    B = 0           #  wave direction [rad]

    # for infinite depth
    g = 9.81            # gravitational constant (m/s^2)
    k = w**2/g          # wave number infinite depth (rad^2/m)
    lam = int(2*np.pi/k)    # wavelength infinite depth (m)

    # defining mesh
    body = cpt.FloatingBody(cpt.meshes.predefined.rectangles.mesh_parallelepiped(size=(th, wi, h), 
                                                                                 resolution=(nt, nw, nh), 
                                                                                 center=(x, y, z), 
                                                                                 name='rect'))
    body.keep_immersed_part()
    body.center_of_mass=(0,0,-h/2)
    body.add_all_rigid_body_dofs()
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()

    array = body + body.translated((xtrans[0],ytrans[0],0),name='2') + body.translated((xtrans[1],ytrans[1],0),name='3')
    array.add_all_rigid_body_dofs()
    #array.show_matplotlib()

    # solving diffraction
    solver = cpt.BEMSolver()
    diff_prob = cpt.DiffractionProblem(body=array, wave_direction=B, omega=w)
    diff_result = solver.solve(diff_prob,keep_details=(True))

    # generating wave elevation dataset
    x1, x2, nx, y1, y2, ny = -res*lam, res*lam, res*lam, -res*lam, res*lam, res*lam
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    fse = solver.compute_free_surface_elevation(grid, diff_result)
    incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)
    total = fse + incoming_fse
    kd = total/incoming_fse
    return kd, total, incoming_fse, lam