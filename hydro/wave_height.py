def wave_height(total, incoming_fse,lam,res,xtrans,ytrans,farm):
    import numpy as np
    from scipy.signal import find_peaks
    import pandas as pd
    import warnings
    warnings.filterwarnings("ignore", category=np.ComplexWarning)
    ##################################################################
    # Extract relevant columns
    numerical_margin = 10/2
    zinc_up = incoming_fse[int((res/2)*lam-1),:int(((res/2)*lam)-numerical_margin)]
    zinc_down = incoming_fse[int((res/2)*lam-1),int(((res/2)*lam)+numerical_margin):]

    z_up = total[int((res/2)*lam-1),:int(((res/2)*lam)-numerical_margin)]
    z_down = total[int((res/2)*lam-1),int(((res/2)*lam)+numerical_margin):]

    avg_H_zup = np.array([np.mean(abs(z_up))])
    avg_H_zincup = np.array([np.mean(abs(zinc_up))])
    avg_H_zdown = np.array([np.mean(abs(z_down))])
    avg_H_zincdown = np.array([np.mean(abs(zinc_down))])

    if farm == True:
        z_upWEC1 = total[int(((res/2)*lam-1) + ytrans[0]/2),:int(((res/2)*lam)+(xtrans[0]/2)-numerical_margin)]
        z_upWEC3 = total[int(((res/2)*lam-1) + ytrans[1]/2),:int(((res/2)*lam)+(xtrans[1]/2)-numerical_margin)]

        z_downWEC1 = total[int(((res/2)*lam-1) + ytrans[0]/2),int(((res/2)*lam) + (xtrans[0]/2)+numerical_margin):]
        z_downWEC3 = total[int(((res/2)*lam-1) + ytrans[1]/2),int(((res/2)*lam) + (xtrans[1]/2)+numerical_margin):]

        avg_H_zup = np.array([np.mean(abs(z_upWEC1)),np.mean(abs(z_up)),np.mean(abs(z_upWEC3))])
        avg_H_zdown = np.array([np.mean(abs(z_downWEC1)),np.mean(abs(z_down)),np.mean(abs(z_downWEC3))])

    ref_H = abs(avg_H_zup/avg_H_zincup) - 1
    trans_H = abs(avg_H_zdown/avg_H_zincdown)
    EB1 = trans_H**2 + ref_H**2     # energy balance
    EB2 = (trans_H + ref_H)**2      # energy balance
    
    return ref_H, trans_H, EB1, EB2