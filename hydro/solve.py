def hydro(array,B,depth,w,farm,controls):
    import capytaine as cpt
    import numpy as np
    import PTO

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
    RAO_vals = abs((RAO.values))            # this is essentially the true pitch amplitude

    if controls == True:
        RAO_controlled = PTO.RAO(diff_prob,diff_result,dataset,array,w,farm)
        RAO_vals = RAO_controlled
    return diff_result,rad_result,RAO_vals,lam

def elevation(res,lam,diff_result,rad_result,RAO_vals,farm,rad,controls,N,attenuator):
    import numpy as np
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
    import capytaine as cpt
    import matplotlib.pyplot as plt

    solver = cpt.BEMSolver()

    # defining the computational grid and preparing post-process data
    x1, x2, nx, y1, y2, ny = -res*lam, res*lam, res*lam, -res*lam, res*lam, res*lam
    grid = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    diffraction = solver.compute_free_surface_elevation(grid, diff_result)  # wave el due to diffraction
    multiplications = []
    if attenuator == True:
        single = 4
        multiple = single*N
    else:
        single = 1
        multiple = 3
    if farm == False:
        for i in range(single):
            if controls == True:
                mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * RAO_vals[i]
                multiplications.append(mult_result)
            else:
                mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * RAO_vals[0,i]
                multiplications.append(mult_result)
    else:
        for i in range(multiple):
            if controls == True:
                mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * RAO_vals[i]
                multiplications.append(mult_result)
            else:
                mult_result = solver.compute_free_surface_elevation(grid, rad_result[i]) * RAO_vals[0,i]
                multiplications.append(mult_result)
    radiation = sum(multiplications)
    if rad == False:
        radiation = 0
    incoming_fse = airy_waves_free_surface_elevation(grid, diff_result)     # incident wave el
    total = diffraction + radiation + incoming_fse                          # total wave el

    # import matplotlib.patheffects as path_effects
    # xtrans = np.array([50,0,50])                        # x translation of bodies if farm
    # ytrans = np.array([50,0,-50])
    # # plots
    # Z = np.real(total)
    # X = grid[0]
    # Y = grid[1]
    # plt.pcolormesh(X, Y, Z)
    # plt.xlabel("x")
    # plt.ylabel("y")
    # colorbar = plt.colorbar()
    # colorbar.set_label(r"Total Wave Elevation, $\eta$")
    # # Add markers with black outline
    # plt.scatter(xtrans, ytrans, marker='o', color='red', s=100, edgecolor='black', linewidth=1.5)

    # # Add arrow with black outline
    # #plt.arrow(-150, 150, 50, 0, color='black', width=0.3, head_width=7, head_length=7)
    # #plt.arrow(-150, 150, 50, 0, color='red', width=0.2, head_width=5, head_length=5)

    # # Add text with black outline
    # #text = plt.text(-100, 130, 'Incident Waves', color='red', fontsize=12, ha='center', va='center')
    # #text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

    # plt.tight_layout()
    # plt.savefig('PA_OS.pdf')
    # plt.show()

    return total, incoming_fse