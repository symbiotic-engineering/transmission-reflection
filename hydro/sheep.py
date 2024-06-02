def wec_run(w,breakwtr,point_absorber,oscillating_surge,attenuator,farm,controls):
    ################ SCRIPT FOR CALCULATING AT COEFFS ############
    # run all Kt and Kr calcs for any body here
    import body         # "body" contains functions for PA, OSWEC, attenuator, and breakwater initialization
    import solve        # this solves the hydrodynamics
    import wave_height  # this function finds Kt and Kr based on wave elevation
    import numpy as np
    import matplotlib.pyplot as plt 

    # setting farm to "False" simulates an isolated body. if it is
    # set to "True", you will simulate three bodies
    if farm == False:
        index = 1
    else:
        index = 3

    B = 0                                           # wave direction [rad]
    depth = 40                                      # average water depth at southfork
    xtrans = np.array([0,0])                        # x translation of bodies if farm
    ytrans = np.array([50,-50])                     # y translation of bodies if farm

    # attenuator properties (because it has to be modeled a little differently)
    N = 3                           # number of attenuators
    D = 30                          # distance btwn cylinders in a single attenuator

    Kr_H = [[] for _ in range(index)]       # initializing reflection coeff
    Kt_H = [[] for _ in range(index)]       # initializing transmission coeff

    w_vals = []                             # for storing omega values

    # here, we change the resolution of the computational grid of
    # the free surfuce wrt to wavelength. this keeps computational time reasonable
    # for long wavelength and accuracy up to par for short wavelengths
    for w in w:
        res = 2.5
        # this where the code generates the body based on which you chose
        if breakwtr == True:
            array, rel_dim = body.breakwater(xtrans,ytrans,farm)
            rad = False               # rad only false for breakwater case
        else:
            rad = True
        if point_absorber == True:
            array, rel_dim = body.PA(xtrans,ytrans,farm)
        if oscillating_surge == True:
            array, rel_dim = body.OSWEC(xtrans,ytrans,farm)
        if attenuator == True:
            array, rel_dim = body.attenuator(xtrans,ytrans,farm,D)

        # this is where the code solves hydrodynamics
        diff_result,rad_result,RAO_vals,lam = solve.hydro(array,B,depth,w,farm,controls)
        # this is where the code solves for wave elevation
        total, incoming_fse, x1, x2, nx, y1, y2, ny = solve.elevation(res,lam,diff_result,rad_result,RAO_vals,farm,rad,controls,N,attenuator,rel_dim)

        # this is where the code calculates your reflection and transmission coefficients
        ref, trans, EB, KD = wave_height.wave_height(total, incoming_fse,lam,xtrans,ytrans,farm,rel_dim,nx,ny,x1,x2,y1,y2)
        print('energy balance',EB)
        print('Kt',trans)
        print('Kr',ref)
        print('dissipation',KD)
        w_vals.append(w)
        for i in range(index):
            Kr_H[i].append(ref[i])
            Kt_H[i].append(trans[i])
        
    return Kt_H, Kr_H, w_vals


# this is where i've been generating my Kt(omega)/Kr(omega) plots
import numpy as np
import matplotlib.pyplot as plt

file_path = '/mnt/c/Users/ov162/transmission-reflection/hydro/figures/'
file_name = 'PA_reg_uncontrolled.pdf'

w = np.array([0.7,0.8,0.9,1.0,1.1,1.2,1.3]) #([0.5,0.65,0.75,0.85,0.95,1.047,1.2,1.25,1.3])   # wave frequency

Kt_H, Kr_H, w_vals = wec_run(w,breakwtr=False,point_absorber=True,oscillating_surge=False,
                             attenuator=False,farm=True,controls=False)

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