def hydro(array,B,depth,w,farm,controls):
    import capytaine as cpt
    import numpy as np
    import controls

    g = 9.81            # gravitational constant (m/s^2)
    k = w**2/g      # wave number infinite depth (rad^2/m)
    lam = int(2*np.pi/k)    # wavelength infinite depth (m)
    
    # solving hydrodynamics
    solver = cpt.BEMSolver()
    diff_prob = cpt.DiffractionProblem(body=array, wave_direction=B, water_depth=depth,omega=w)
    diff_result = solver.solve(diff_prob,keep_details=(True))
    rad_prob = [
        cpt.RadiationProblem(body=array, radiating_dof=dof, water_depth=depth,omega=w)
        for dof in array.dofs
        ]
    rad_result = solver.solve_all(rad_prob,keep_details=(True))
    dataset = cpt.assemble_dataset(rad_result + [diff_result])

    RAO = cpt.post_pro.rao(dataset, wave_direction=B, dissipation=None, stiffness=None)
    RAO_vals = abs(np.array((RAO.values)))            # this is essentially the true pitch amplitude

    RAO_controlled = controls.RAO(diff_prob,diff_result,dataset,array,w,farm)
    if controls == True:
        RAO_vals = RAO_controlled
    return diff_result,rad_result,RAO_vals,lam

def elevation(res,lam,diff_result,rad_result,RAO_vals,farm,rad):
    import numpy as np
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
    import capytaine as cpt
    solver = cpt.BEMSolver()

    # defining the computational grid and preparing post-process data
    x1, x2, nx, y1, y2, ny = -res*lam, res*lam, res*lam, -res*lam, res*lam, res*lam
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    diffraction = solver.compute_free_surface_elevation(grid, diff_result)  # wave el due to diffraction
    multiplications = []
    if farm == False:
        for i in range(1):
            mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * RAO_vals[0,i]
            multiplications.append(mult_result)
    else:
        for i in range(3):
            mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * RAO_vals[0,i]
            multiplications.append(mult_result)
    radiation = sum(multiplications)
    if rad == False:
        radiation = 0
    incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)     # incident wave el
    total = diffraction + radiation + incoming_fse                          # total wave el
    kd = total/incoming_fse

    return kd, total, incoming_fse