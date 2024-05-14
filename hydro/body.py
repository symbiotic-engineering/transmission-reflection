def PA(xtrans,ytrans,farm):
    # define body
    import capytaine as cpt

    # initializing parameters
    r = 10     # radius [m]
    l = 5
    x = 0
    y = 0
    nr = 20         # number of panels along radius
    ntheta = 10     # number of panels in theta direction
    nz = 5         # number of panels in z-direction
    z = 0

    # defining mesh
    body = cpt.FloatingBody(mesh=cpt.mesh_vertical_cylinder(length=l, radius=r,center=(x,y,z),
                                                            resolution=(nr,ntheta,nz),name='cyl'))
    body.keep_immersed_part()
    body.center_of_mass=(0,0,z)
    body.add_all_rigid_body_dofs()
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    body.keep_only_dofs(dofs='Heave')

    # create array
    array = body + body.translated((xtrans[0],ytrans[0],0),name='2') + body.translated((xtrans[1],ytrans[1],0),name='3')
    array.add_all_rigid_body_dofs()
    array.keep_only_dofs(dofs=['cyl__Heave','2__Heave','3__Heave'])

    if farm == False:
        array = body
    
    return array

def OSWEC(xtrans,ytrans,farm):
    import capytaine as cpt
    # initializing parameters
    wi = 25         # width of flap [m]
    th = 1          # thickness of flap [m]
    h = 19          # height of flap [m], draft = 16 m
    x = 0
    y = 0
    draft = 16 # m
    z = 0.5*h-draft # box center [m]
    cog = -11.4     # center of gravity [m]
    nw = 10         # number of panels along width (x)
    nt = 3          # number of panels along thickness (y)
    nh = 10         # number of panels along height (z)

    # defining mesh
    body = cpt.FloatingBody(cpt.meshes.predefined.rectangles.mesh_parallelepiped(size=(th, wi, h), 
                                                                                resolution=(nt, nw, nh), 
                                                                                center=(x, y, z),
                                                                                name='rect'))
    body.keep_immersed_part()
    body.center_of_mass=(0,0,cog)
    body.add_all_rigid_body_dofs()
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    body.keep_only_dofs(dofs='Pitch')

    # create array
    array = body + body.translated((xtrans[0],ytrans[0],0),name='2') + body.translated((xtrans[1],ytrans[1],0),name='3')
    array.add_all_rigid_body_dofs()
    array.keep_only_dofs(dofs=['rect__Pitch','2__Pitch','3__Pitch'])
    
    if farm == False:
        array = body
    
    return array

def breakwater(xtrans,ytrans,farm):
    import capytaine as cpt

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

    if farm == False:
        array = body
    # array.show_matplotlib()
    return array