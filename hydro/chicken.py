def wec_run(w):
    ################ FUNCTION FOR USE WHEN RUNNING ALL WITH SWAN ############
    # run all Kt and Kr calcs for any body here
    # "body" contains functions for PA, OSWEC, and breakwater initialization
    import body
    import attenuators
    import PA_OSWEC
    import wave_height
    import kd_post
    import solve
    import numpy as np
    import matplotlib.pyplot as plt

    # if you want an array, farm = True
    # for single body, farm = False
    farm = False
    rad = False                                      # only false for breakwater case

    B = 0                                           # wave direction [rad]
    depth = 40                                      # average water depth at southfork
    xtrans = np.array([0,0])
    ytrans = np.array([50,-50])

    ########attenuator ONLY (right now)###############
    N = 3                           # number of bodies
    D = 30                          # distance btwn bodies
    ##################################################

    Kr_H = []
    Kt_H = []
    Kr_K = []
    Kt_K = []
    w_vals = []
    for w in w:
        if w < 0.8:
            res = 2
        else: 
            res = 2                                   # resolution factor of grid wrt lambda
        array = body.breakwater(xtrans,ytrans,farm)
        kd, total, incoming_fse, lam = solve.solver(array,B,depth,w,res,farm,rad)
        #kd, total, incoming_fse, lam = attenuators.lpf(w,res,N,D,farm)
        ref_H, trans_H, EB1, EB2 = wave_height.wave_height(total, incoming_fse,lam, res)
        ref_K, trans_K, EB1, EB2 = kd_post.disturbance(kd, lam, res)
        Kr_H.append(ref_H)
        Kt_H.append(trans_H)
        Kr_K.append(ref_K)
        Kt_K.append(trans_K)
        w_vals.append(w)

    return Kr_H, Kt_H, Kr_K, Kt_K