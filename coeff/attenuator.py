###MAKING MAJOR ASSUMPTION THAT BODY HAS ONE PITCH DOF AND NOT INDIVIDUAL FOR
###EACH JOINT. WILL TRY AND CHANGE LATER ON####
def bodysolver(w):

    import capytaine as cpt
    from scipy.linalg import block_diag
    import numpy as np
    from scipy.interpolate import griddata

    # initializing parameters
    r = 1.75     # radius [m]
    l = 29
    x = -50
    y = 0
    nr = 10        # number of panels along radius
    ntheta = 10     # number of panels in theta direction
    nz = 10         # number of panels in z-direction
    z = 0
    B = 0      # wave direction [rad]
    depth = 40      # average water depth at southfork

    # defining mesh
    body = cpt.FloatingBody(mesh=cpt.mesh_horizontal_cylinder(length=l, radius=r,center=(x,y,z),
                                                              resolution=(nr,ntheta,nz),name='cyl'))
    body.keep_immersed_part()
    body.center_of_mass=(0,0,0)
    body.add_all_rigid_body_dofs()
    body.inertia_matrix = body.compute_rigid_body_inertia()
    body.hydrostatic_stiffness = body.compute_hydrostatic_stiffness()
    print(body.dofs.keys())
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
    print(pitch_RAO)
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
    print('e',e)

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
    # multiplications = []
    # for i in range(4):
    #     mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * pitch_RAO[0,i]
    #     multiplications.append(mult_result)
    # radiation = sum(multiplications)
    return radiation, lam, grid