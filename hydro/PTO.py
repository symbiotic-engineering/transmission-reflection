'''This file extracts the hydro coefficients, froude-krylov force, and 
diffraction force from the diffraction problem and uses them to compute:
1. the exciting force
2. the optimal damping and stiffness PTO terms (to resonate)
3. the controlled response amplitude operator (RAO)'''

def RAO(diff_prob,diff_result,dataset,array,w,farm,char_dim,point_absorber,reactive):

    from capytaine.bem.airy_waves import froude_krylov_force
    import numpy as np

    # extract hydro coeffs
    # to include off-diagonals
    A = np.squeeze(np.array([[dataset['added_mass'].sel(radiating_dof=effecting, influenced_dof=effected) for effecting in array.dofs] for effected in array.dofs]))
    B = np.squeeze(np.array([[dataset['radiation_damping'].sel(radiating_dof=effecting, influenced_dof=effected) for effecting in array.dofs] for effected in array.dofs]))
    K = array.hydrostatic_stiffness.values
    M = array.inertia_matrix.values

    # extract forces and compute exciting force
    FK = np.array([froude_krylov_force(diff_prob)[dof] for dof in array.dofs])
    dif = np.array([diff_result.forces[dof] for dof in array.dofs])
    ex_force = FK + dif
    
    # Define simple optimal PTO damping and stiffness
    # for reactive control:
    if reactive:
        B_pto = B
        K_pto = w**2*(M+A)-K  
    else:
        # for damping only:
        B_pto = (B**2 + ( w*(M+A) - (K/w) )**2)**0.5
        K_pto = 0 

    # FOR INCLUDING OFF-DIAGONALS
    inertia = M + A 
    resistance = B + B_pto
    reactance = K + K_pto
    H = -(w**2)*inertia - 1j*w*resistance + reactance 

    if farm:
        RAO_controlled = np.linalg.solve(H,ex_force).ravel()
    else:
        RAO_controlled = ex_force/H
    print('RAO_controlled',RAO_controlled)

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
    if farm:
        power = 0.5*np.diag(B_pto)*(abs(RAO_controlled*w*1j))**2        # [kW]
    else:
        power = 0.5*B_pto*(abs(RAO_controlled*w*1j))**2        # [kW]

    # power available in wave
    rho = 1025                  # [kg/m^3] density of sea water
    g = 9.81                    # [m/s^2] gravitational constant
    power_avail = (rho * g**2 * amplitude**2) / (4 * w)      # [kW/m]

    CW = power/power_avail               # [m]
    CWR = CW / char_dim                  # unitless
    print('capture width ratio',CWR)

    return RAO_controlled, CWR
