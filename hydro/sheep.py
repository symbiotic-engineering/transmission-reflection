def wec_run(w,breakwtr,point_absorber,oscillating_surge,attenuator,farm,controls):
    ################ SCRIPT FOR CALCULATING AT COEFFS ############
    # run all Kt and Kr calcs for any body here
    # "body" contains functions for PA, OSWEC, attenuator, and breakwater initialization
    import body
    import wave_height
    import solve
    import numpy as np
    import matplotlib.pyplot as plt 

    if farm == False:
        index = 1
    else:
        index = 3

    B = 0                                           # wave direction [rad]
    depth = 40                                      # average water depth at southfork
    xtrans = np.array([50,50])                        # x translation of bodies if farm
    ytrans = np.array([50,-50])                     # y translation of bodies if farm

    ######## attenuator properties ###############
    N = 3                           # number of attenuators
    D = 30                          # distance btwn cylinders in a single attenuator
    ##################################################

    Kr_H = [[] for _ in range(index)]
    Kt_H = [[] for _ in range(index)]

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
        else:
            rad = True
        if point_absorber == True:
            array = body.PA(xtrans,ytrans,farm)
        if oscillating_surge == True:
            array = body.OSWEC(xtrans,ytrans,farm)
        if attenuator == True:
            array = body.attenuator(xtrans,ytrans,farm,D)

        # this is where the code solves hydrodynamics
        diff_result,rad_result,RAO_vals,lam = solve.hydro(array,B,depth,w,farm,controls)
        # this is where the code solves for wave elevation
        total, incoming_fse = solve.elevation(res,lam,diff_result,rad_result,RAO_vals,farm,rad,controls,N,attenuator)

        # this is where the code calculates your reflection and transmission coefficients
        ref_H, trans_H, EB1, EB2 = wave_height.wave_height(total, incoming_fse,lam, res,xtrans,ytrans,farm)
        w_vals.append(w)

        for i in range(index):
            Kr_H[i].append(ref_H[i])
            Kt_H[i].append(trans_H[i])
        

    return Kt_H, Kr_H, w_vals


import numpy as np
import matplotlib.pyplot as plt

file_path = '/mnt/c/Users/ov162/transmission-reflection/hydro/figures/'
file_name = 'test.pdf'

w = np.array([0.5,0.65,0.75,0.85,0.95,1.047,1.2,1.25,1.3])   # wave frequency

Kt_H, Kr_H, w_vals = wec_run(w,breakwtr=False,point_absorber=True,oscillating_surge=False,
                             attenuator=False,farm=False,controls=False)

cud_colors = ['#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7', '#000000', '#8B4513']
linestyles = ['-', '--', ':', '-.', '-', '--', ':', '-.']

for i, kt_h_values in enumerate(Kt_H):
    plt.plot(w_vals, kt_h_values, marker='o', label=f'$K_t$ for body {i+1}', 
             color=cud_colors[i % len(cud_colors)], linestyle=linestyles[i % len(linestyles)])
for i, kr_h_values in enumerate(Kr_H):
    plt.plot(w_vals, kr_h_values, marker='x', label=f'$K_r$ for body {i+1}', 
             color=cud_colors[i+1 % len(cud_colors)], linestyle=linestyles[i % len(linestyles)])

plt.legend()
plt.xlabel('$\omega$ [rad/s]')
plt.ylabel('Coefficient Value')
plt.savefig(f'{file_path}{file_name}')
plt.show()