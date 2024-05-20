def wec_run(w,breakwtr,point_absorber,oscillating_surge,attenuator,farm):
    ################ SCRIPT FOR CALCULATING AT COEFFS ############
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

    if farm == False:
        index = 1
    else:
        index = 3

    B = 0                                           # wave direction [rad]
    depth = 40                                      # average water depth at southfork
    xtrans = np.array([0,0])                        # x translation of bodies if farm
    ytrans = np.array([50,-50])                     # y translation of bodies if farm

    ######## attenuator ONLY (right now) ###############
    N = 3                           # number of bodies
    D = 30                          # distance btwn bodies
    ##################################################

    Kr_H = [[] for _ in range(index)]
    Kt_H = [[] for _ in range(index)]
    # Kr_K = [[] for _ in range(index)]
    # Kt_K = [[] for _ in range(index)]

    w_vals = []                             # for storing omega values
    for w in w:
        if w < 0.8:
            res = 2
        else: 
            res = 3                                   # resolution factor of grid wrt lambda

        # this where the code generates the body based on which you chose
        if breakwtr == True:
            array = body.breakwater(xtrans,ytrans,farm)
            rad = False               # rad only false for breakwater case
            kd, total, incoming_fse, lam = solve.solver(array,B,depth,w,res,farm,rad)
        else:
            rad = True
        if point_absorber == True:
            array = body.PA(xtrans,ytrans,farm)
            kd, total, incoming_fse, lam = solve.solver(array,B,depth,w,res,farm,rad)
        if oscillating_surge == True:
            array = body.OSWEC(xtrans,ytrans,farm)
            kd, total, incoming_fse, lam = solve.solver(array,B,depth,w,res,farm,rad)
        if attenuator == True:
            kd, total, incoming_fse, lam = attenuators.lpf(w,res,N,D,farm)

        # this is where the code calculates your reflection and transmission coefficients
        ref_H, trans_H, EB1, EB2 = wave_height.wave_height(total, incoming_fse,lam, res,xtrans,ytrans,farm)
        # ref_K, trans_K, EB1, EB2 = kd_post.disturbance(kd, lam, res)
        w_vals.append(w)

        for i in range(index):
            Kr_H[i].append(ref_H[i])
            Kt_H[i].append(trans_H[i])
    return Kt_H, Kr_H, w_vals


import numpy as np
import matplotlib.pyplot as plt

file_path = '/mnt/c/Users/ov162/transmission-reflection/hydro/figures/'
file_name = 'atten_array_coeff.pdf'

w = np.array([0.5,0.65,0.75,0.85,0.95,1.047,1.2])   # wave frequency

Kt_H, Kr_H, w_vals = wec_run(w,breakwtr=False,point_absorber=False,oscillating_surge=False,attenuator=True,farm=True)

cud_colors = ['#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7', '#000000', '#8B4513']
linestyles = ['-', '--', ':', '-.', '-', '--', ':', '-.']

for i, kt_h_values in enumerate(Kt_H):
    plt.plot(w_vals, kt_h_values, marker='o', label=f'$K_t$ for body {i+1}', 
             color=cud_colors[i % len(cud_colors)], linestyle=linestyles[i % len(linestyles)])
for i, kr_h_values in enumerate(Kr_H):
    plt.plot(w_vals, kr_h_values, marker='x', label=f'$K_r$ for body {i+1}', 
             color=cud_colors[i % len(cud_colors)], linestyle=linestyles[i % len(linestyles)])

plt.legend()
plt.xlabel('$\omega$ [rad/s]')
plt.ylabel('Coefficient Value')
plt.savefig(f'{file_path}{file_name}')
plt.show()