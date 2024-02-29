def bodysolver(w):

    import capytaine as cpt
    from scipy.linalg import block_diag
    import numpy as np
    from scipy.interpolate import griddata

    # initializing parameters
    wi = 25         # width of flap [m]
    th = 1          # thickness of flap [m]
    h = 19          # height of flap [m], draft = 16 m
    x = 0
    y = 0
    draft = 16 # m
    z = 0.5*h-draft # box center [m]
    cog = -11.4     # center of gravity [m]
    nw = 15         # number of panels along width (x)
    nt = 3          # number of panels along thickness (y)
    nh = 15         # number of panels along height (z)
    B = 0           # wave direction [rad]
    depth = 40      # average water depth at southfork

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

    # solving hydrodynamics
    solver = cpt.BEMSolver()
    diff_prob = cpt.DiffractionProblem(body=body, wave_direction=B, water_depth=depth,omega=w)
    diff_result = solver.solve(diff_prob,keep_details=(True))
    rad_prob = cpt.RadiationProblem(body=body, radiating_dof='Pitch', water_depth=depth,omega=w)
    rad_result = solver.solve(rad_prob,keep_details=(True))
    dataset = cpt.assemble_dataset([rad_result] + [diff_result])
    RAO = cpt.post_pro.rao(dataset, wave_direction=B, dissipation=None, stiffness=None)
    pitch_RAO = np.array(np.abs(RAO.values))            # this is essentially the true pitch amplitude

    return dataset, rad_result, rad_prob, diff_prob, diff_result, body, pitch_RAO
    
def exc_force(diff_prob, rad_result, body, w):
    import numpy as np
    from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity
    
    rho = 1020      # density of sea water
    # Read mesh properties
    faces_centers = body.mesh.faces_centers
    faces_normals = body.mesh.faces_normals
    faces_areas = body.mesh.faces_areas

    # Get potentials
    phi_inc = airy_waves_potential(faces_centers, diff_prob)
    v_inc = airy_waves_velocity(faces_centers, diff_prob)
    phi_rad = rad_result.potential                            # radiation potential on body

    # Indirect computation from the radiation solution, via the Haskind relation
    integrand = - (phi_inc * faces_normals[:,2]
                   - phi_rad * np.diag(v_inc@faces_normals.T))
    e = 1j * w * rho * np.sum(integrand*faces_areas)               # force from the Haskind relation

    return e

def hydro(dataset, body):
    # hydro coeffs
    damp = (dataset['radiation_damping'].sel(radiating_dof=['Pitch'],
                                            influenced_dof=['Pitch'])).values
    add = (dataset['added_mass'].sel(radiating_dof=['Pitch'],
                                    influenced_dof=['Pitch'])).values
    add = add[0][0][0]
    damp = damp[0][0][0]
    stiff = body.hydrostatic_stiffness.values
    mass = body.inertia_matrix.values
    return damp, add, stiff, mass

def phi_rad(w, rad_result, diff_result, pitch_RAO):
    import capytaine as cpt
    import numpy as np
    g = 9.81            # gravitational constant (m/s^2)
    lam = int(2*np.pi/(w**2/g))   # wavelength infinite depth (m)
    x1, x2, nx, y1, y2, ny = -4*lam, 4*lam, 4*lam, -4*lam, 4*lam, 4*lam
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    solver = cpt.BEMSolver()
    radiation = (solver.compute_free_surface_elevation(grid, rad_result))*pitch_RAO  # wave el due to pitch radiation
    diffraction = solver.compute_free_surface_elevation(grid,diff_result)
    tot = radiation + diffraction       # not including incident
    return tot, lam, grid