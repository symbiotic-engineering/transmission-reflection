def RAO(diff_prob,diff_result,dataset,array,w,farm):
    # simple controls model to obtain controlled RAO for radiation elevation

    from capytaine.bem.airy_waves import airy_waves_potential, airy_waves_velocity, froude_krylov_force
    import numpy as np

    B = np.array([dataset['radiation_damping'].sel(radiating_dof=dofs,influenced_dof=dofs) for dofs in array.dofs])
    A = np.array([dataset['added_mass'].sel(radiating_dof=dofs, influenced_dof=dofs) for dofs in array.dofs])
    K = array.hydrostatic_stiffness.values
    M = array.inertia_matrix.values

    # heave exciting forces
    FK = np.array([froude_krylov_force(diff_prob)[dof] for dof in array.dofs])
    dif = np.array([diff_result.forces[dof] for dof in array.dofs])
    ex_force = FK + dif
    
    # Define simple PTO damping and stiffness
    B_pto = B
    K_pto = w**2*(M+A)-K  
    # B_pto = (B**2+(w*(M+A)-K/w)**2)**0.5
    # K_pto = 0 
    
    # WEC motion (complex) 
    RAO_controlled = abs((np.diag(ex_force/((-w**2)*(M+A) - (B + B_pto)*w*1j + K + K_pto)))) 
         
    return RAO_controlled

