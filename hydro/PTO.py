'''This file extracts the hydro coefficients, froude-krylov force, and 
diffraction force from the diffraction problem and uses them to compute:
1. the exciting force
2. the optimal damping and stiffness PTO terms (to resonate)
3. the controlled response amplitude operator (RAO)'''

def RAO(diff_prob,diff_result,dataset,array,w,farm):

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
    B_pto = B
    K_pto = w**2*(M+A)-K  

    # for damping only:
    # B_pto = (B**2+(w*(M+A)-K/w)**2)**0.5
    # K_pto = 0 
    
    # WEC motion (complex) 
    # RAO is used in further computation
    RAO_controlled = (np.diag(ex_force/((-w**2)*(M+A) - (B + B_pto)*w*1j + K + K_pto)))
    print('rao controlled',RAO_controlled)

    # power is computed but not used
    power = 0.5*np.diag(B_pto)*(abs(RAO_controlled))**2*w**2
    print('power',power)
         
    return RAO_controlled

