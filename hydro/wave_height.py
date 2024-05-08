def wave_height(total, incoming_fse,lam,res):
    import numpy as np
    from scipy.signal import find_peaks
    import pandas as pd
    import warnings
    warnings.filterwarnings("ignore", category=np.ComplexWarning)
    ##################################################################
    # Extract relevant columns
    z_up = total[int((res/2)*lam-1),:int((res/2)*lam)]
    zinc_up = incoming_fse[int((res/2)*lam-1),:int((res/2)*lam)]
    z_down = total[int((res/2)*lam-1),int((res/2)*lam):]
    zinc_down = incoming_fse[int((res/2)*lam-1),int((res/2)*lam):]

    avg_H_zup = np.mean(abs(z_up))
    avg_H_zincup = np.mean(abs(zinc_up))
    avg_H_zdown = np.mean(abs(z_down))
    avg_H_zincdown = np.mean(abs(zinc_down))

    ref_H = abs(avg_H_zup/avg_H_zincup) - 1
    trans_H = abs(avg_H_zdown/avg_H_zincdown)
    EB1 = trans_H**2 + ref_H**2     # energy balnce
    EB2 = (trans_H + ref_H)**2      # energy balance
    
    return ref_H, trans_H, EB1, EB2