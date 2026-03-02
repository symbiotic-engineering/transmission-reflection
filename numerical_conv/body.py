'''for mesh convergence use only
'''

def initialize(xtrans,ytrans,w,x_center,point_absorber,nr,ntheta,nz,nt,nw,nh):
    import capytaine as cpt
    import matplotlib.pyplot as plt
    import numpy as np
    import logging
    logging.getLogger('capytaine').setLevel(logging.ERROR)
    
    rho_w = 1000  
    
    ##### initializing PA parameters #####
    r,l = 14.6/2, 4.445                                 # radius [m], length [m] #10.5, 6          
    m_PA = 431000
    COB = m_PA/(2*np.pi*r**2*rho_w)                     # center of buoyancy based on mass and radius
    Pdraft = 2*COB
    z = 0.5 * l - Pdraft                                # body center position           
    PA_cog = (3.23469) - Pdraft                         # center of gravity (from top of rack)
    PA_com = np.array([0,0,PA_cog])                     # center of mass           
    #nr, ntheta, nz = 16, 32, 12                         # panels in each direction 

    rel_dim = r                           # dimension relevant for computing Kt and Kr while avoiding body location (orthogonal to wave)
    A_W = np.pi*(r**2)                    # maximum horizontal cross-sectional area of device [m^2]
    char_dim = np.sqrt((4*A_W)/np.pi)     # dimension relevant for computing power available in wave taken from Babarit CW classification

    ##### initializing OSWEC parameters #####
    wi,th,h = 18, 1.905, 10.8                     # width, thickness, and height of flap [m]
    draft = 9.2                                   # [m]
    z_flap = 0.5 * h - draft                      # box center [m]
    OS_cog = -7.6                                 # center of gravity [m] from aisha's calculations
    OS_com = np.array([0,0,OS_cog])
    #nw,nt,nh = 30, 10, 18                         # number of panels along width (x), thickness (y), and height (z)

    rel_dim = th/2                                # dimension relevant for computing Kt and Kr while avoiding body location (orthogonal to wave)
    char_dim = wi                                 # dimension relevant for computing power available in wave


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
    x_start, y_start = x_center, 0                          
    Dx, Dy = xtrans, ytrans

    if point_absorber:
        PA_positions = (x_start, y_start, z)
        PA_names = 'cyl1'
        PA_bodies =  create_floating_body(cpt.mesh_vertical_cylinder, {**pa_mesh_args, 'center': PA_positions}, PA_com, 'Heave', PA_names) 

        array = PA_bodies
        # print(PA_bodies)
        # array.show_matplotlib()
        # plt.savefig('PAmesh.pdf')
        # plt.clf
    else:
        OS_positions = (x_start, y_start, z_flap)
        OS_names = 'rect'
        OS_bodies = create_floating_body(cpt.meshes.predefined.rectangles.mesh_parallelepiped, {**oswec_mesh_args, 'center': OS_positions}, OS_com, 'Surge', OS_names)
        array = OS_bodies
        array.show_matplotlib()
        plt.savefig('OSmesh.pdf')
        plt.clf

    return array, rel_dim, char_dim