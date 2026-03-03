'''MESH CONVERGENCE VERSION'''

def RAO(diff_prob,diff_result,dataset,array,w,char_dim,point_absorber,controls):

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
    
    if controls:
        B_pto = (B**2 + ( w*(M+A) - (K/w) )**2)**0.5
        K_pto = 0 
    else:
        B_pto = 0
        K_pto = 0

    # FOR INCLUDING OFF-DIAGONALS
    inertia = M + A 

    resistance = B + B_pto
    reactance = K + K_pto
    H = -(w**2)*inertia - 1j*w*resistance + reactance 

    RAO = ex_force/H

    # power produced by WEC, used to find CWR
    power = 0.5*(B_pto)*(abs(RAO*w*1j))**2        # [W]

    # power available in wave
    rho = 1025                  # [kg/m^3] density of sea water
    g = 9.81                    # [m/s^2] gravitational constant
    amplitude = 1
    power_avail = (rho * g**2 * amplitude**2) / (4 * w)      # [kW/m]

    CW = power/power_avail               # [m]
    CWR = CW / char_dim                  # unitless

    return RAO, CWR