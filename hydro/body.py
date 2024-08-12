'''this file contains initialization for the following bodies:
1. Point Absorber (based off of Reference Model 3 from NREL)
2. Oscillating Surge WEC (based off Reference Model 5 from NREL)
3. Floating Breakwater (based off Naval Standard)
4. Attenuator (based off Pelamis device)

This script will generate a single body and an array. It will
fix the bodies such that the only degree of freedom considered is the
one in which they extract power. Mesh sizes were determined from the
mesh_conv.py script. The "rel_dim" is the "relevant dimension," used
in the wave_height.py script to ensure the Kt and Kr coefficients 
are not calculated over a space the body occupies.'''

def PA(xtrans,ytrans,farm,w,x_center):
    import capytaine as cpt
    import matplotlib.pyplot as plt
    import numpy as np
    import logging
    logging.getLogger('capytaine').setLevel(logging.ERROR)

    # initializing parameters
    r,l = 10.5, 6              # radius [m], length (5) [m]
    x, y, z = x_center, 0, 0        # body center position           
    cog = -0.3*(l/2)         # center of mass as reported by WECSim
    nr, ntheta, nz = 14, 25, 9      # panels in each direction 19, 6

    # dimension relevant for computing Kt and Kr while avoiding
    # body location (orthogonal to wave)
    rel_dim = r

    # dimension relevant for computing power available in wave
    # taken from Babarit CW classification
    A_W = np.pi*(r**2)                # maximum horizontal cross-sectional area of device [m^2]
    char_dim = np.sqrt((4*A_W)/np.pi) # ends up just being the diameter [m]

    # budal's power limit (heaveing, axisymmetric)
    rho = 1025                  # density of sea water [kg/m^3]
    g = 9.81                    # gravitational constant [m/s^2]
    A = 1                       # unit wave amplitude [m]
    V = np.pi * r**2 * l        # volume of device [m^3]
    budal_limit = (1/4)*rho*g*A*w*V     # Budal's power limit for heaving axisymmetric body [kg-m^2/s^3]

    # defining the mesh and create floating body
    cylinder_mesh = cpt.mesh_vertical_cylinder(length=l, radius=r,center=(x,y,z),
                                                            resolution=(nr,ntheta,nz),
                                                            name='cyl')
    #lid_mesh = cylinder_mesh.generate_lid(z=-0.1)
    body = cpt.FloatingBody(mesh=cylinder_mesh)#,lid_mesh=lid_mesh)
    body.keep_immersed_part()               # clips body for computation
    body.center_of_mass=np.array([0,0,cog])             # defines center of mass (required)
    body.add_all_rigid_body_dofs()          # capytaine requirement
    body.inertia_matrix = body.compute_rigid_body_inertia()             # compute inertia matrix (required)
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()   # compute hydrostatic stiffness (required)
    body.keep_only_dofs(dofs='Heave')
    #body.show_matplotlib()

    # create array
    array = body + body.translated((xtrans[0],ytrans[0],0),name='2') + body.translated((xtrans[1],ytrans[1],0),name='3')
    array.add_all_rigid_body_dofs()
    array.keep_only_dofs(dofs=['cyl__Heave','2__Heave','3__Heave'])     # this fixes the body in all other DOFs

    if farm == False:
        array = body
        #array.show_matplotlib()
        #plt.savefig('pa_panels.pdf')

    return array, rel_dim, char_dim, budal_limit

def OSWEC(xtrans, ytrans, farm, w,x_center):
    import capytaine as cpt
    import matplotlib.pyplot as plt
    import numpy as np
    import logging
    logging.getLogger('capytaine').setLevel(logging.ERROR)

    # initializing parameters
    wi, th, h = 18, 1, 16         # width, thickness, and height of flap [m]
    draft = h - 3                    # draft [m]
    x, y, z = x_center, 0, 0.5*h-draft   # postion of body center
    cog = -0.7125*draft                   # center of gravity [m] (71.25% of the draft)
    nt, nh, nw = 2, 24, 27        # number of panels in each direction     

    # dimension relevant for computing Kt and Kr while avoiding
    # body location (orthogonal to wave)
    rel_dim = th/2 

    # dimension relevant for computing power available in wave
    char_dim = wi

    # budal's power limit (pitching, non-axisymmetric)
    rho = 1025                  # density of sea water [kg/m^3]
    A = 1                       # unit amplitude [m]
    budal_limit = (1/24) * np.pi * rho * A * w**3 * wi**3 * th

    # defining mesh
    body = cpt.FloatingBody(cpt.meshes.predefined.rectangles.mesh_parallelepiped(size=(th, wi, h), 
                                                                                resolution=(nt, nw, nh), 
                                                                                center=(x, y, z),
                                                                                name='rect'))
    body.keep_immersed_part()
    body.center_of_mass=np.array([0,0,cog])
    body.add_all_rigid_body_dofs()
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    body.keep_only_dofs(dofs='Pitch')
    #body.show_matplotlib()

    # create array
    array = body + body.translated((xtrans[0],ytrans[0],0),name='2') + body.translated((xtrans[1],ytrans[1],0),name='3')
    array.add_all_rigid_body_dofs()
    array.keep_only_dofs(dofs=['rect__Pitch','2__Pitch','3__Pitch'])
    
    if farm == False:
        array = body
    return array, rel_dim, char_dim, budal_limit

def breakwater(xtrans,ytrans,farm,x_center):
    import capytaine as cpt
    import matplotlib.pyplot as plt
    import numpy as np
    import logging
    logging.getLogger('capytaine').setLevel(logging.ERROR)
    # initializing parameters
    wi, th, h = 20, 5, 2         # width, thickness, and height of box [m]
    x, y, z = x_center, 0, -0.5         # box center (x = -25 for staggered array, x = 0 otherwise))
    cog = -0.162                 # 10.8% of the draft, equivalent to PA cog
    nw, nt, nh = 12, 6, 5        # number of panels in each direction

    # dimension relevant for computing Kt and Kr while avoiding
    # body location (orthogonal to wave)
    rel_dim = th/2
    # dimension relevant for computing power available in wave
    char_dim = wi

    # defining mesh
    body = cpt.FloatingBody(cpt.meshes.predefined.rectangles.mesh_parallelepiped(size=(th, wi, h), 
                                                                                 resolution=(nt, nw, nh), 
                                                                                 center=(x, y, z), 
                                                                                 name='rect'))
    body.keep_immersed_part()
    body.center_of_mass=np.array([0,0,cog])
    body.add_all_rigid_body_dofs()
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    #body.show_matplotlib()

    array = body + body.translated((xtrans[0],ytrans[0],0),name='2') + body.translated((xtrans[1],ytrans[1],0),name='3')
    array.add_all_rigid_body_dofs()
    #array.show_matplotlib()
    #plt.savefig('break_arraystag.pdf')

    if farm == False:
        array = body

    return array, rel_dim, char_dim

def attenuator(xtrans,ytrans,farm,w,x_center):
    import capytaine as cpt
    import matplotlib.pyplot as plt
    import numpy as np
    import logging
    logging.getLogger('capytaine').setLevel(logging.ERROR)

    # initializing parameters
    r, l = 2, 14                     # radius, length (of one cylinder) [m]
    total_length = l*2 + 1
    x = -(1/2)*(1 + l)                # formula for two body
    x, y, z = x + x_center, 0, 0                 # body center of first cylinder
    nr, ntheta, nz = 4, 17, 16        # number of panels in each direction
    cog = -(r/2)*0.108                # 10.8% of the draft, equivalent to PA cog, tough to find for pelamis
    D = l + 1                         # distance btwn cylinder centers (in the single attenuator)
    
    # dimension relevant for computing Kt and Kr while avoiding
    # body location (orthogonal to wave)
    rel_dim = int(total_length/2)

    # dimension relevant for computing power available in wave
    char_dim = l

    # budal's power limit (pitching, axisymmetric)
    rho = 1025                                  # density of sea water [kg/m^3]
    A = 1                                       # unit wave amplitude [m]
    budal_limit = (rho*A/8) * np.pi**2 * r**2 * (l**2 + r**2) * w**3

    # defining mesh
    body = cpt.FloatingBody(mesh=cpt.mesh_horizontal_cylinder(length=l, radius=r,center=(x,y,z),
                                                            resolution=(nr,ntheta,nz),name='cyl'))
    body.keep_immersed_part()
    body.center_of_mass=np.array([0,0,cog])
    body.add_all_rigid_body_dofs()
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()

    array = body + body.translated((D,0,0),name='1b') + body.translated((xtrans[0],ytrans[0],0),name='2a') + body.translated((D + xtrans[0],ytrans[0],0),name='2b') + body.translated((xtrans[1],ytrans[1],0),name='3a') + body.translated((D+xtrans[1],ytrans[1],0),name='3b')
    array.keep_only_dofs(dofs=['cyl__Pitch','1b__Pitch','2a__Pitch','2b__Pitch','3a__Pitch','3b__Pitch'])


    if farm == False:
        array = body + body.translated((D,0,0),name='1b')
        array.keep_only_dofs(dofs=['cyl__Pitch','1b__Pitch'])
        #array.show_matplotlib()
        #plt.savefig('atten_body.pdf')
    return array, rel_dim, char_dim, budal_limit