def disturbance(kd,lam,res):
    import pandas as pd
    import numpy as np

    xR = np.linspace(-1*int(res*lam),0,int((res/2)*lam))
    xD = np.linspace(0,int(res*lam),int((res/2)*lam))

    #########################################################################################################

    # Extract relevant columns
    kd_up = kd[int((res/2)*lam)-1,:int((res/2)*lam)]
    avg_kd_up = np.mean(abs(kd_up))
    ref_K = abs(avg_kd_up) - 1

    ######################################## FOR DOWNSTREAM ###########################################
    kd_down = kd[int((res/2)*lam)-1,int((res/2)*lam):]
    avg_distance_kd_down = np.mean(abs(kd_down))
    trans_K= abs(avg_distance_kd_down)

    EB1 = trans_K**2 + (ref_K - 1)**2       # energy balance check
    EB2 = (trans_K + (ref_K - 1))**2       # energy balance check
    return ref_K, trans_K, EB1, EB2