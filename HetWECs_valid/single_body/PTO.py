''' This script is adapted from the TEAMER HetWECs project by VITALE.
Edits made for validation purposes by VITALE 07/10/2025
'''
def RAO(diff_prob,diff_result,dataset,body,w,reactive,b,PA,B_difference):
    # optimal controls model to obtain controlled RAO for radiation elevation based on Falnes theory
    from capytaine.bem.airy_waves import froude_krylov_force
    import numpy as np
    import capytaine as cpt
    import friction
    
    # extract hydro coeffs
    A = np.squeeze(np.array([[dataset['added_mass'].sel(radiating_dof=effecting, influenced_dof=effected) for effecting in body.dofs] for effected in body.dofs]))
    B = np.squeeze(np.array([[dataset['radiation_damping'].sel(radiating_dof=effecting, influenced_dof=effected) for effecting in body.dofs] for effected in body.dofs]))
    K = body.hydrostatic_stiffness.values
    M = body.inertia_matrix.values

    # correcting capytaine hydro coeffs
    if PA == False:
        K = 1.153               # [kg-m^2/s^2]
        M = 0.0405              # [kg-m^2]
        #A = A*2.95              # gets natural period of 3.0426 s
    else:
        M = 3.448
        A = A   #*1.0725

    nat_per = (2*np.pi)/(np.sqrt(K/(A+M)))
    #print('nat per',nat_per)
    # print('input per',np.pi*2/w)
    # print('diff', (np.pi*2/w) - nat_per)
    # if np.abs((np.pi*2/w) - nat_per) < 0.0001:
    #     print('YAAAAAAAAAAHOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO')

    # extract forces and compute exciting force
    FK = np.array([froude_krylov_force(diff_prob)[dof] for dof in body.dofs])
    dif = np.array([diff_result.forces[dof] for dof in body.dofs])
    ex_force = FK + dif

    #surge_force = FK[1] + dif[1]
    
    # Define simple optimal PTO damping and stiffness for reactive control:
    if reactive:
        B_pto = B
        #K_pto = w**2*(M+A)-K  
        #B_pto = 0
        K_pto = 0
    else:
        B_pto = 0
        K_pto = 0
    
    inertia = M + A 
    resistance = B + B_pto + B_difference
    reactance = K + K_pto
    H = -(w**2)*inertia - 1j*w*resistance + reactance 
    #RAO_controlled = np.linalg.solve(H,ex_force).ravel()

    RAO_controlled = ex_force/H
    k = w**2/9.81
    power = (0.5*abs(B_pto)*abs((w*H*0.03)**2))
    print('power',np.abs(power))
    
    # loop for testing if budal limit is violated and for adding viscous drag to flap 
    if PA:
        B_viscous = 0
        # # using coloumb friction term from Mi et al. 2024 (almost exactly the same as the one in my
        # # 'friction.py' script in the archive folder)
        # T_f = 3*w**2                              # initial guess [N-m]
        # B_friction = (4 * T_f) / (np.pi * w * np.abs(RAO_controlled))

        # # B_friction = friction.slidefriction(B,M,K,A,w)

        amplitude = 1.00
        body_velocity = np.abs(RAO_controlled) * 1j * w
        if np.abs(body_velocity) > amplitude * w:
            print('Budal limit violated at w = ',w)

    else:
        B_viscous = 0
        # gam_nl = 0.00
        # B_viscous = (((8 * np.abs(RAO_controlled) * w) / (3 * np.pi)) * gam_nl)
        # B_friction = 0
        # resistance = B + B_pto + B_viscous + B_friction + B_difference
        # H = -(w**2)*inertia - 1j*w*resistance + reactance 
        # RAO_controlled = ex_force/H


    return RAO_controlled, ex_force, A, B, B_viscous, nat_per


    # # using coloumb friction term from Mi et al. 2024 (almost exactly the same as the one in my
    # # 'friction.py' script in the archive folder)
    # T_f = 0.75*(ex_force)                               # initial guess [N-m]
    # B_friction = (4 * T_f) / (np.pi * w * np.abs(RAO_controlled))
    # resistance = B + B_pto + B_viscous + B_friction
    # H = -(w**2)*inertia - 1j*w*resistance + reactance 
    # RAO_controlled = ex_force/H

    # # adding reynolds's viscous friction
    # mu = 0.03                               # guess for viscous friction coeff in bearing
    # F_vf = np.abs(ex_force*mu*(RAO_controlled* 1j * w))
    # print('visc fric bearing',(F_vf))
    # RAO_controlled = (ex_force - F_vf)/H