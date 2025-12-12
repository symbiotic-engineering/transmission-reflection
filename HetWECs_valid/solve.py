def hydro(array,B,depth,w,reactive,B_diffPA,B_diffOS,megRAOexp,joRAOexp,virRAOexp,franRAOexp,
            Bd_vir,Bd_fran,Bd_meg,Bd_jo,PA,wg_amp_exp,het):
    import capytaine as cpt
    import numpy as np
    import PTO
    
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

    RAO, ex_force, added_mass, damping, B_PTOemp, mech_power, dissipative_power = PTO.RAO(diff_prob,diff_result,dataset,
                                                                                    array,w,reactive,B,B_diffPA,B_diffOS,
                                                                                    megRAOexp,joRAOexp,virRAOexp,franRAOexp,
                                                                                    Bd_vir,Bd_fran,Bd_meg,Bd_jo,PA,wg_amp_exp,het)

    return RAO, diff_result, rad_result, ex_force, added_mass, damping, B_PTOemp, mech_power, dissipative_power

def elevation(res,diff_result,rad_result,RAO_vals,xlocs,ylocs,wg_amp_exp):
    import numpy as np
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
    import capytaine as cpt
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata

    solver = cpt.BEMSolver()

    # defining the computational grid and preparing post-process data
    x1, x2, y1, y2 = 12 - (13.086 + 1.55), 18 - (13.086 + 1.55), -1.5 - (-0.1410 - 0.50), 1.5 - (-0.1410 - 0.50) #5, 25, -2.5, 4
    ny = int(res*2*(abs(y1)+y2))
    nx = int(res*(abs(x1)+x2))
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    diffraction = solver.compute_free_surface_elevation(grid, diff_result) # wave elevation due to diffraction

    multiplications = []

    for i in range(4):
        mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * RAO_vals[i]
        multiplications.append(mult_result)
    radiation = sum(multiplications)
    incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)                  # incident wave elevation
    total = (incoming_fse + diffraction + radiation)* wg_amp_exp                          # total wave elevation(diffraction + incoming_fse + radiation)

    # Interpolate 'total' onto the gauge points
    X = grid[0]
    Y = grid[1]
    points = np.column_stack((xlocs, ylocs))
    elevation_at_gauges = griddata((X.ravel(), Y.ravel()), total.ravel(), points, method='linear')
    # total_at_gauges now contains the values of 'total' at the gauge locations specified by gauge_x and gauge_y

    return total, incoming_fse, grid, radiation, diffraction,elevation_at_gauges