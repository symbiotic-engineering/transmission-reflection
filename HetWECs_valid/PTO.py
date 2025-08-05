''' This script is adapted from the TEAMER HetWECs project by VITALE.
Edits made for validation purposes by VITALE 07/10/2025
'''
def RAO(diff_prob,diff_result,dataset,array,w,reactive,b):
    # optimal controls model to obtain controlled RAO for radiation elevation based on Falnes theory

    from capytaine.bem.airy_waves import froude_krylov_force
    import numpy as np
    import capytaine as cpt
    
    # extract hydro coeffs
    # to include off-diagonals
    A = np.squeeze(np.array([[dataset['added_mass'].sel(radiating_dof=effecting, influenced_dof=effected) for effecting in array.dofs] for effected in array.dofs]))
    B = np.squeeze(np.array([[dataset['radiation_damping'].sel(radiating_dof=effecting, influenced_dof=effected) for effecting in array.dofs] for effected in array.dofs]))
    K = array.hydrostatic_stiffness.values
    M = array.inertia_matrix.values

    # for i, dof in enumerate(array.dofs):
    #     if dof.startswith(('rect', 'rect2')):
    #         K = 1.82       # [kg-m^2/s^2]
    #         M = 0.0461     # [kg-m^2]

    # extract forces and compute exciting force
    FK = np.array([froude_krylov_force(diff_prob)[dof] for dof in array.dofs])
    dif = np.array([diff_result.forces[dof] for dof in array.dofs])
    ex_force = dif #FK + dif
    
    # Define simple optimal PTO damping and stiffness for reactive control:
    if reactive:
        B_pto = B
        K_pto = w**2*(M+A)-K  
    else:
        B_pto = 0
        K_pto = 0

    # FOR INCLUDING OFF-DIAGONALS
    inertia = M + A 
    resistance = B + B_pto
    reactance = K + K_pto
    H = -(w**2)*inertia - 1j*w*resistance + reactance 
    RAO_controlled = np.linalg.solve(H,ex_force).ravel()

    # B_original = B
    
    # for i, dof in enumerate(array.dofs):
    #     if dof.startswith(('rect', 'rect2')):

    #         # accounting for viscous damping coeff for flap
    #         tolerance = 1e-3  # Define your tolerance level
    #         max_iterations = 100000  # Maximum iterations to prevent infinite loop
    #         i = 0  # Starting index for RAO_controlled

    #         gam_nl = 0.1  # Initial nonlinear viscous damping coeff guess
    #         previous_RAO_controlled = RAO_controlled  # To store RAO from the previous iteration
    #         iteration = 0

    #         while iteration < max_iterations:
    #             # Calculate B_viscous, resistance, reactance, H, and RAO_controlled
    #             B_viscous = ((8 * RAO_controlled * w) / (3 * np.pi)) * gam_nl
    #             B = B_original + B_viscous  # Use the original B for each iteration
    #             resistance = B + B_pto
    #             H = -(w ** 2) * inertia - 1j * w * resistance + reactance
    #             RAO_controlled = np.linalg.solve(H, ex_force).ravel()
    #             cond_H = np.linalg.cond(H)
    #             if cond_H > 1e10:
    #                 print(f"Warning: Ill-conditioned H with cond={cond_H}")

    #             # Check for convergence
    #             if previous_RAO_controlled is not None:
    #                 if np.all(abs(RAO_controlled - previous_RAO_controlled) < tolerance):
    #                     print(f"Convergence achieved with gam_nl = {gam_nl}")
    #                     break

    #             previous_RAO_controlled = RAO_controlled
    #             gam_nl += 0.0001  # You can adjust the increment step based on the behavior
    #             iteration += 1

    #         if iteration == max_iterations:
    #             print("Max iterations reached without convergence.")
    # print('post loop',RAO_controlled)

    # amplitude = 0.03  # unit wave amplitude [m]
 
    # for i in range(np.size(RAO_controlled)):
    #     if (i==2) or (i==3):
    #         print('running loop')
    #         body_velocity = RAO_controlled[i] * 1j * w
            
    #         while np.abs(body_velocity) > amplitude * w:

    #             # Update RAO_controlled
    #             RAO_controlled[i] = 0.99 * (np.real(RAO_controlled[i]) + 1j * np.imag(RAO_controlled[i]))
                
    #             # Recalculate body velocity
    #             body_velocity = RAO_controlled[i] * 1j * w

    #         print('in loop',RAO_controlled)
    #   print('post loop',RAO_controlled)


    return RAO_controlled, ex_force