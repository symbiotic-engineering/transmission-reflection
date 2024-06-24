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

    for w in w:

        # this is where you set the grid resolution based on wavelength.
        # shorter wavelengths require finer mesh resolution based on 
        # mesh convergence study.
        if w < 0.95:
            res = 1
        elif w == 1.0:
            res = 2.0
        elif w == 1.1:
            res = 2.5
        elif w == 1.25:
            res = 2.5
        elif w == 1.3:
            res = 3.5

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
import csv

# file_path = '/mnt/c/Users/ov162/transmission-reflection/hydro/figures/'
# file_name = 'PA_reg_uncontrolled.pdf'

w = np.array([0.7,0.8,0.9,1.0,1.1,1.25,1.3])   # wave frequency

Kt_H, Kr_H, w_vals = wec_run(w,breakwtr=False,point_absorber=False,oscillating_surge=False,
                             attenuator=True,farm=False,controls=True)

print('frequency',w_vals)
print('transmission coeff',Kt_H)
print('reflection coeff',Kr_H)

cud_colors = ['#009E73','#D55E00','#E69F00', '#56B4E9','#0072B2', '#CC79A7', '#000000','#F0E442']
linestyles = ['-', '--', ':', '-.', '-', '--', ':', '-.']

for i, kt_h_values in enumerate(Kt_H):
    plt.plot(w_vals, kt_h_values, marker='s', markersize='10',label=f'$K_t$ body {i+1}',
             color=cud_colors[i % len(cud_colors)], linestyle=linestyles[i % len(linestyles)],linewidth=3)
for i, kr_h_values in enumerate(Kr_H):
    plt.plot(w_vals, kr_h_values, marker='d', markersize='10', label=f'$K_r$ body {i+1}', 
             color=cud_colors[i+1 % len(cud_colors)], linestyle=linestyles[i % len(linestyles)],linewidth=3)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('$\omega$ [rad/s]', fontsize=20)
plt.ylabel('Coefficient Value', fontsize=20)
plt.tight_layout()
print('beep')
plt.savefig('atten_damp.pdf')
# plt.savefig(f'{file_path}{file_name}')
print('boop')
#plt.show()

# this is to save the data if you just want to fiddle with plots
file_path = 'atten_damp.csv'
with open(file_path, mode='w', newline='') as file:
    writer = csv.writer(file)
    header = ['Omega'] + [f'Kt_H_{i+1}' for i in range(len(Kt_H))] + [f'Kr_H_{i+1}' for i in range(len(Kr_H))]
    writer.writerow(header)
    for i in range(len(w_vals)):
        row = [w_vals[i]] + [Kt_H[j][i] for j in range(len(Kt_H))] + [Kr_H[j][i] for j in range(len(Kr_H))]
        writer.writerow(row)