''' this function is adapted from the TEAMER HetWECs repository. the devices modeled in this script
correspond to physical prototypes that were tested in the O.H. Hinsdale Driectional Wave Basin.
Validation of the model will be conducted against experimental data gathered from that campaign.
The physical devices were built at a 1:50 scale, so the meshes will be scaled as well.
'''
def initialize():
    # hydro
    import capytaine as cpt
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import csv
    import random

    # Froude scale factor used to scale prototypes
    scale = 50
    Dx = Dy = 40/scale

    # point absorber dimensions
    r,l = (14.6/2)/scale, 4.445/scale       # radius [m], length [m]
    rho_w = 1000            
    m_PA = 3.448
    COB = m_PA/(2*np.pi*r**2*rho_w)         # center of buoyancy based on mass and radius
    Pdraft = 2*COB
    z = 0.5 * l - Pdraft                    # body center position           
    nr, ntheta, nz = 16, 22, 12             # panels in each direction
    PA_cog = (3.23469/scale) - Pdraft       # center of gravity (from top of rack)
    PA_com = np.array([0,0,PA_cog])         # center of mass

    # OSWEC dimensions
    wi,th,h = 18/scale, 1.905/scale, 10.8/scale   # width, thickness, and height of flap [m]
    draft = 0.184                                 # m
    z_flap = 0.5 * h - draft                      # box center [m]
    OS_cog = -0.152                               # center of gravity [m] from aisha's calculations
    OS_com = np.array([0,0,OS_cog])
    nw,nt,nh = 30, 10, 18                         # number of panels along width (x), thickness (y), and height (z)
    ####### ACTUAL OSWEC MASS = 1.625 KG

    def create_floating_body(mesh_func, mesh_args, center_of_mass, dofs, name):
        body = cpt.FloatingBody(mesh=mesh_func(**mesh_args, name=name))     # create floating body
        body.keep_immersed_part()                                           # trim mesh 
        body.center_of_mass = center_of_mass                                # define center of mass
        body.add_all_rigid_body_dofs()                                      # generate degrees of freedom
        body.inertia_matrix = body.compute_rigid_body_inertia()             # compute inertia matrix
        body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()   # compute hydrostatic stiffness
        body.keep_only_dofs(dofs=dofs)                                      # choosing relevant dofs (defined in function inputs)
        return body

    # Define mesh arguments for point absorbers and oscillating surges
    pa_mesh_args = {
        'length': l,
        'radius': r,
        'resolution': (nr, ntheta, nz)
    }

    oswec_mesh_args = {
        'size': (th, wi, h),
        'resolution': (nt, nw, nh)
    }

    # Create point absorbers and oscillating surges
    OS_positions = [(14.9900 - 0.25, -0.1410 - 0.25, z_flap), (14.9900 - 0.25, -0.1410+Dy - 0.25, z_flap)]
    PA_positions = [(14.9900 + Dx - 0.25, -0.1410 + (1/2)*Dy - 0.25, z), (14.9900 + Dx - 0.25, -0.1410 + (3/2)*Dy - 0.25, z)]
    OS_names = ['rect', 'rect2']
    PA_names = ['cyl3', 'cyl4']

    # generate meshed bodies
    PA_bodies = [
        create_floating_body(cpt.mesh_vertical_cylinder, {**pa_mesh_args, 'center': pos}, PA_com, 'Heave', name) 
        for pos, name in zip(PA_positions, PA_names)
    ]

    OS_bodies = [
        create_floating_body(cpt.meshes.predefined.rectangles.mesh_parallelepiped, {**oswec_mesh_args, 'center': pos}, OS_com, 'Surge', name)
        for pos, name in zip(OS_positions, OS_names)
    ]

    # Assign PA and OS variables for het4b configuration
    OS, OS2 = OS_bodies
    PA3, PA4 = PA_bodies

    # create meshed array
    array = OS + OS2 + PA3 + PA4
    # array.show_matplotlib()
    # plt.savefig('array.pdf')
    # plt.clf


    return array