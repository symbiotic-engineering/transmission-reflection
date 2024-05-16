################ SCRIPT FOR ONLY LOOKING AT COEFFS ############
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
farm = True
rad = True               # rad only false for breakwater case

w = np.array([0.5,0.65,0.75,0.85,0.95,1.047,1.2])   # wave frequency
B = 0                                           # wave direction [rad]
depth = 40                                      # average water depth at southfork
xtrans = np.array([0,0])
ytrans = np.array([50,-50])

########attenuator ONLY (right now)###############
N = 3                           # number of bodies
D = 30                          # distance btwn bodies
##################################################
if farm == False:
    Kr_H = []
    Kt_H = []
    Kr_K = []
    Kt_K = []
else:
    Kr_H = [[] for _ in range(3)]
    Kt_H = [[] for _ in range(3)]
    Kr_K = [[] for _ in range(3)]
    Kt_K = [[] for _ in range(3)]

w_vals = []
for w in w:
    if w < 0.8:
        res = 2
    else: 
        res = 3                                   # resolution factor of grid wrt lambda
    #array = body.OSWEC(xtrans,ytrans,farm)
    #kd, total, incoming_fse, lam = solve.solver(array,B,depth,w,res,farm,rad)
    kd, total, incoming_fse, lam = attenuators.lpf(w,res,N,D,farm)
    ref_H, trans_H, EB1, EB2 = wave_height.wave_height(total, incoming_fse,lam, res,xtrans,ytrans,farm)
    ref_K, trans_K, EB1, EB2 = kd_post.disturbance(kd, lam, res)
    if farm == False:
        Kr_H.append(ref_H)
        Kt_H.append(trans_H)
        Kr_K.append(ref_K)
        Kt_K.append(trans_K)
        w_vals.append(w)
    else:
        for i in range(3):
            Kr_H[i].append(ref_H[i])
            Kt_H[i].append(trans_H[i])
            #Kr_K[i].append(ref_K[i])
            #Kt_K[i].append(trans_K[i])
        w_vals.append(w)

if farm==False:
    plt.plot(w_vals, Kt_H, color='blue', marker='o', label='$K_t$ from H')
    plt.plot(w_vals, Kr_H, color='red', marker='o', label='$K_r$ from H')
    plt.plot(w_vals, Kt_K, color='blue', marker='o', linestyle=':', label='$K_t$ from $K_d$')
    plt.plot(w_vals, Kr_K, color='red', marker='o', linestyle=':', label='$K_r$ from $K_d$')
    plt.legend()
    plt.xlabel('$\omega$ [rad/s]')
    plt.ylabel('Coefficient Value')
    plt.savefig('/mnt/c/Users/ov162/transmission-reflection/hydro/figures/atten_coeffs.pdf')
    plt.show()
else:
    linestyles = ['-', '--', ':', '-.', '+', 'x', '*']
    for i, kt_h_values in enumerate(Kt_H):
        plt.plot(w_vals, kt_h_values, marker='o', label=f'$K_t$ for body {i+1}', 
                 linestyle=linestyles[i % len(linestyles)])
    for i, kr_h_values in enumerate(Kr_H):
        plt.plot(w_vals, kr_h_values, marker='o', label=f'$K_r$ for body {i+1}', 
                 linestyle=linestyles[i % len(linestyles)])

    plt.legend()
    plt.xlabel('$\omega$ [rad/s]')
    plt.ylabel('Coefficient Value')
    plt.savefig('/mnt/c/Users/ov162/transmission-reflection/hydro/figures/atten_array_coeff.pdf')
    plt.show()
