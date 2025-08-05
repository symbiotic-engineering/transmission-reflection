def hydro(body,B,depth,w,reactive,PA):
    import capytaine as cpt
    import numpy as np
    import PTO
    
    # solving hydrodynamics
    solver = cpt.BEMSolver()
    diff_prob = cpt.DiffractionProblem(body=body, wave_direction=B, water_depth=depth,omega=w)
    diff_result = solver.solve(diff_prob,keep_details=(True))
    rad_prob = [
        cpt.RadiationProblem(body=body, radiating_dof=dof, water_depth=depth,omega=w)
        for dof in body.dofs
        ]
    rad_result = solver.solve_all(rad_prob,keep_details=(True))
    dataset = cpt.assemble_dataset(rad_result + [diff_result])

    RAO, ex_force = PTO.RAO(diff_prob,diff_result,dataset,body,w,reactive,B,PA)

    return RAO, diff_result, rad_result, ex_force

def elevation(res,diff_result,rad_result,RAO_vals):
    import numpy as np
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
    import capytaine as cpt
    import matplotlib.pyplot as plt

    solver = cpt.BEMSolver()

    # defining the computational grid and preparing post-process data
    x1, x2, y1, y2 = -200, 200, -75, 125
    ny = int(res*(abs(y1)+y2))
    nx = int(res*(abs(x1)+x2))
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    diffraction = solver.compute_free_surface_elevation(grid, diff_result)  # wave elevation due to diffraction

    multiplications = []

    for i in range(4):
        mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * RAO_vals[i]
        multiplications.append(mult_result)
    radiation = sum(multiplications)
    incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)     # incident wave elevation
    total = diffraction + incoming_fse + radiation                          # total wave elevation

    return total, incoming_fse, grid, radiation, diffraction