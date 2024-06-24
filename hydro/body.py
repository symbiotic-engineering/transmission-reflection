def PA(xtrans,ytrans,farm):
    import capytaine as cpt
    import matplotlib.pyplot as plt
    import numpy as np
    # initializing parameters
    r = 10              # radius [m]
    l = 5               # length [m]
    x = 0               # x-position of body center
    y = 0               # y-position of body center
    z = 0               # z-position of body center
    cog = -0.72         # center of mass as reported by WECSim
    nr = 12
    ntheta = 19          # number of panels in theta direction (was 10)
    nz = 4              # number of panels in z-direction (was 5)
    rel_dim = r + abs(xtrans[0])


    # defining the mesh and create floating body
    body = cpt.FloatingBody(mesh=cpt.mesh_vertical_cylinder(length=l, radius=r,center=(x,y,z),
                                                            resolution=(nr,ntheta,nz),name='cyl'))
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
        rel_dim = r
        array = body
    return array, rel_dim

def OSWEC(xtrans, ytrans, farm):
    import capytaine as cpt
    import matplotlib.pyplot as plt
    import numpy as np

    # initializing parameters
    wi = 25         # width of flap [m]
    th = 1          # thickness of flap [m]
    h = 19          # height of flap [m]
    draft = 16      # draft [m]
    x = 0           # x-postion of body center
    y = 0           # y-position of body center
    z = 0.5*h-draft # z-position of body center [m]
    cog = -11.4     # center of gravity [m]
    nt = 3         
    nh = 10              
    nw = 18
    rel_dim = th + abs(xtrans[0])

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
        rel_dim = th
    return array, rel_dim

def breakwater(xtrans,ytrans,farm):
    import capytaine as cpt
    import matplotlib.pyplot as plt
    import numpy as np
    # initializing parameters
    wi = 20         # width of box [m]
    th = 5          # thickness of box [m]
    h = 2           # height of box [m]
    x = 0
    y = 0
    z = -0.5        # box center [m]
    cog = -0.162    # 10.8% of the draft, equivalent to PA cog
    nw = 12         # number of panels along width (x)
    nt = 6         # number of panels along thickness (y)
    nh = 5          # number of panels along height (z)
    rel_dim = th + abs(xtrans[0])

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

    if farm == False:
        array = body
        rel_dim = th
    return array, rel_dim

def attenuator(xtrans,ytrans,farm,D):
    # D = distance btwn cylinders (in the single attenuator) = 30
    import capytaine as cpt
    import matplotlib.pyplot as plt
    import numpy as np
    # initializing parameters
    r = 1.75     # radius [m]
    l = 29
    x = -50
    y = 0
    nr = 1        # number of panels along radius (was 1)
    ntheta = 15    # number of panels in theta direction (was 6)
    nz = 41         # number of panels in length-direction (41)
    z = 0
    cog = -0.189   # 10.8% of the draft, equivalent to PA cog, tough to find for pelamis
    rel_dim = int(((l*4)/2) + abs(xtrans[0]))

    # defining mesh
    body = cpt.FloatingBody(mesh=cpt.mesh_horizontal_cylinder(length=l, radius=r,center=(x,y,z),
                                                            resolution=(nr,ntheta,nz),name='cyl'))
    body.keep_immersed_part()
    body.center_of_mass=np.array([0,0,cog])
    body.add_all_rigid_body_dofs()
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()

    array = body + body.translated((D,0,0),name='1b') + body.translated((D*2,0,0),name='1c') + body.translated((D*3,0,0),name='1d') + body.translated((xtrans[0],ytrans[0],0),name='2a') + body.translated((D + xtrans[0],ytrans[0],0),name='2b') + body.translated((D*2+xtrans[0],ytrans[0],0),name='2c') + body.translated((D*3+xtrans[0],ytrans[0],0),name='2d')+ body.translated((xtrans[1],ytrans[1],0),name='3a') + body.translated((D+xtrans[1],ytrans[1],0),name='3b') + body.translated((D*2+xtrans[1],ytrans[1],0),name='3c') + body.translated((D*3+xtrans[1],ytrans[1],0),name='3d')
    array.keep_only_dofs(dofs=['cyl__Pitch','1b__Pitch','1c__Pitch','1d__Pitch','2a__Pitch','2b__Pitch',
                                '2c__Pitch','2d__Pitch','3a__Pitch','3b__Pitch','3c__Pitch','3d__Pitch'])


    if farm == False:
        array = body + body.translated((D,0,0),name='1b') + body.translated((D*2,0,0),name='1c') + body.translated((D*3,0,0),name='1d')
        array.keep_only_dofs(dofs=['cyl__Pitch','1b__Pitch','1c__Pitch','1d__Pitch'])
        rel_dim = ((l*4)/2)
        #array.show_matplotlib()
    return array, rel_dim