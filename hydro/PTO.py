'''This file extracts the hydro coefficients, froude-krylov force, and 
diffraction force from the diffraction problem and uses them to compute:
1. the exciting force
2. the optimal damping and stiffness PTO terms (to resonate)
3. the controlled response amplitude operator (RAO)'''

def RAO(diff_prob,diff_result,dataset,array,w,farm,char_dim,point_absorber):

    from capytaine.bem.airy_waves import froude_krylov_force
    import numpy as np

    # extract hydro coeffs
    B = np.array([dataset['radiation_damping'].sel(radiating_dof=dofs,influenced_dof=dofs) for dofs in array.dofs])
    A = np.array([dataset['added_mass'].sel(radiating_dof=dofs, influenced_dof=dofs) for dofs in array.dofs])
    K = array.hydrostatic_stiffness.values
    M = array.inertia_matrix.values

    # extract forces and compute exciting force
    FK = np.array([froude_krylov_force(diff_prob)[dof] for dof in array.dofs])
    dif = np.array([diff_result.forces[dof] for dof in array.dofs])
    ex_force = FK + dif
    
    # Define simple optimal PTO damping and stiffness
    # for reactive control:
    #B_pto = B
    #K_pto = w**2*(M+A)-K  

    # for damping only:
    B_pto = (B**2 + ( w*(M+A) - (K/w) )**2)**0.5
    K_pto = 0 
    
    # WEC motion (complex) 
    RAO_controlled = (np.diag(ex_force/((-1*w**2)*(M+A) - (B + B_pto)*w*1j + K + K_pto)))

    amplitude = 1.000  # unit wave amplitude [m]
    if point_absorber:
        body_velocity = RAO_controlled * 1j * w
        RAO_controlled = RAO_controlled.copy()
        # Create a mask for elements where body velocity exceeds the limit
        mask = np.abs(body_velocity) > amplitude * w
        
        while np.any(np.abs(body_velocity[mask]) > amplitude * w):
            # Update RAO_controlled for those specific elements
            RAO_controlled[mask] = 0.95 * (np.real(RAO_controlled[mask]) + 1j * np.imag(RAO_controlled[mask]))
            # Recalculate body velocity for those specific elements
            body_velocity[mask] = RAO_controlled[mask] * 1j * w
            # Update the mask for the next iteration
            mask = np.abs(body_velocity) > amplitude * w

    # power produced by WEC, used to find CWR
    power = 0.5*np.diag(B_pto)*(abs(RAO_controlled*w*1j))**2        # [kW]

    # power available in wave
    rho = 1025                  # [kg/m^3] density of sea water
    g = 9.81                    # [m/s^2] gravitational constant
    power_avail = (rho * g**2 * amplitude**2) / (4 * w)      # [kW/m]

    CW = power/power_avail               # [m]
    CWR = CW / char_dim                  # unitless
    print('capture width ratio',CWR)

    return RAO_controlled, CWR

    # while np.any(power > budal_limit):
    #     # Create a mask for elements that exceed the budal_limit
    #     mask = power > budal_limit
    #     RAO_controlled = RAO_controlled.copy()
    #     # Update RAO_controlled for those specific elements
    #     RAO_controlled = 0.99 * RAO_controlled
    #     # Recalculate power for those specific elements
    #     new_power = 0.5 * B_pto * (abs(RAO_controlled * w * 1j))**2
    #     # Use the mask to update only the relevant elements
    #     power = np.where(mask, new_power, power)

    # if point_absorber:
    #     eps = 1             # factor for heave
    # else:
    #     eps = 2             # factor for pitch and surge

    # k = w**2/g                      # wave number infinite depth (rad^2/m)
    # lam = int(2*np.pi/k)            # wavelength infinite depth (m)
    # CW_max = eps*(lam/(2*np.pi))    # theoretical maximum CW (from Babarit)
    # #print('maximum CW', CW_max/char_dim)

    # while np.any(CW > CW_max):
    #     # Create a mask for elements that exceed the budal_limit
    #     mask = CW > CW_max
    #     RAO_controlled = RAO_controlled.copy()
    #     # Update RAO_controlled for those specific elements
    #     RAO_controlled = 0.99 * RAO_controlled
    #     # Recalculate power for those specific elements
    #     new_power = 0.5 * B_pto * (abs(RAO_controlled * w * 1j))**2
        
    #     # Use the mask to update only the relevant elements
    #     power = np.where(mask, new_power, power)
        
    #     # Update CW and CWR for those specific elements
    #     CW = np.where(mask, power / power_avail, CW)
    #     CWR = np.where(mask, CW / char_dim, CWR)
    
    # # experimental best fit for expected CWR (from Babarit), none given for attenuator
    # if point_absorber:
    #     CWR_exp = (1.3*char_dim + 5.6)/100
    #     CI_95 = 21*np.sqrt(1.1 + ((char_dim-15)**2 / 1090))/100
    #     CWR_exp = CWR_exp - CI_95/1.25       # tuning CWR to obtain EB (/4 for reactive, /1.25 for damp)

    # else:
    #     # for OSWEC, Babarit found char dim did not matter as much
    #     # using similar formula for attenuator bc it pitches, and there's
    #     # no available data otherwise
    #     CWR_exp = 8.5/100
    #     CI_95 = 12*np.sqrt(1.1+ ((char_dim-15)**2 / 1090))/100
    #     CWR_exp = CWR_exp + CI_95       # tuning CWR to obtain EB (OS: /4 for reactive, /2 for damp)
    #                                     # atten: 1 for reactive

    # while np.any(CWR > CWR_exp):
    #     mask = CWR > CWR_exp
    #     RAO_controlled = RAO_controlled.copy()
    #     RAO_controlled = 0.99 * RAO_controlled
    #     new_power = 0.5 * B_pto * (abs(RAO_controlled * w * 1j))**2
    #     power = np.where(mask, new_power, power)
    #     # Update CW and CWR for those specific elements
    #     CW = np.where(mask, power / power_avail, CW)
    #     CWR = np.where(mask, CW / char_dim, CWR)

    #print('final capture width ratio',CWR)
    #print('final saturated power', power)

