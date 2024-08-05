'''This is the run file to run all the near field functions. It runs:
1. body.name() which initializes the body
2. solve.hydro() which solves the body hydrodynamics (and includes the
power take-off function nested within)
3. solve.elevation() which solves for the free surface elevation
4. wave_height.wave_height() which solves for wave height ratios i.e.
the transmission and reflection coefficients.'''

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
    depth = 500                                     # keep deep water assumption for EB
    xtrans = np.array([0,0])                        # x translation of bodies if farm
    ytrans = np.array([50,-50])                     # y translation of bodies if farm

    # attenuator properties (because it has to be modeled a little differently)
    N = 3                           # number of attenuators
    D = 30                          # distance btwn cylinders in a single attenuator

    Kr_H = [[] for _ in range(index)]       # initializing reflection coeff
    Kt_H = [[] for _ in range(index)]       # initializing transmission coeff
    power = [[] for _ in range(index)]      # initializing absorberd power [W/m]

    w_vals = []                             # for storing omega values

    for w in w:
        # this is where you set the grid resolution based on wavelength.
        # shorter wavelengths require finer mesh resolution based on 
        # mesh convergence study.
        # if w < 0.95:
        #     res = 1
        # elif w == 1.0:
        #     res = 2.0
        # elif w == 1.1:
        #     res = 3.5
        # elif w == 1.25:
        #     res = 3.5
        # elif w == 1.3:
        #     res = 3.5
        # else:
        #     res = 2
        if w < 1.0:
            res = 2.0
        else:
            res = 3.5

        # this where the code generates the body based on which you chose
        if breakwtr:
            array, rel_dim, char_dim = body.breakwater(xtrans,ytrans,farm)
            rad = False               # rad only false for breakwater case
        else:
            rad = True
        if point_absorber:
            array, rel_dim, char_dim, budal_limit = body.PA(xtrans,ytrans,farm,w)
        if oscillating_surge:
            array, rel_dim, char_dim, budal_limit = body.OSWEC(xtrans,ytrans,farm,w)
        if attenuator:
            array, rel_dim, char_dim, budal_limit = body.attenuator(xtrans,ytrans,farm,D,w)

        # this is where the code solves hydrodynamics
        diff_result,rad_result,RAO_vals,lam,CWR = solve.hydro(array,B,depth,w,char_dim,farm,controls,point_absorber)

        # this is where the code solves for wave elevation
        total,incoming_fse,x1,x2,nx,y1,y2,ny = solve.elevation(res,lam,diff_result,rad_result,RAO_vals,farm,rad,controls,N,attenuator,rel_dim)

        # this is where the code calculates your reflection and transmission coefficients
        ref,trans,EB,KD,power_abs = wave_height.wave_height(total, incoming_fse,xtrans,ytrans,farm,rel_dim,w,nx,ny,x1,x2,y1,y2)
        if controls == False:
            power_abs = power_abs*0
        print('Kt',trans)
        print('Kr',ref)
        print('dissipation',KD)

        w_vals.append(w)
        for i in range(index):
            Kr_H[i].append(ref[i])
            Kt_H[i].append(trans[i])
            power[i].append(power_abs[i])
        
    return Kt_H, Kr_H, w_vals, power

### COMMENT THIS BIT OUT if you are going to run the big_run.py script ###
# this is where i've been generating my Kt(omega)/Kr(omega) datasets
# for post processing 
'''if you are running the attenuator, you will likely run into the issue of 
overlapping panels in Capytaine. Making slight tweaks to the grid resolution
can mediate this issue.'''

import numpy as np
import matplotlib.pyplot as plt
import csv

w = np.array([0.7,0.8,0.9,1.0,1.1,1.25,1.3])   # wave frequency

# manually set False and True statements according to your goals
# if you set breakwtr=True, you must set controls=False
Kt_H, Kr_H, w_vals, power = wec_run(w,breakwtr=False,point_absorber=True,oscillating_surge=False,
                             attenuator=False,farm=True,controls=False)

# this is to save the data to a .csv file
file_path = 'PA_reg_uncont.csv'
with open(file_path, mode='w', newline='') as file:
    writer = csv.writer(file)
    header = ['Omega'] + [f'Kt_H_{i+1}' for i in range(len(Kt_H))] + [f'Kr_H_{i+1}' for i in range(len(Kr_H))] + [f'power_{i+1}' for i in range(len(power))]
    writer.writerow(header)
    for i in range(len(w_vals)):
        row = [w_vals[i]] + [Kt_H[j][i] for j in range(len(Kt_H))] + [Kr_H[j][i] for j in range(len(Kr_H))] + [power[j][i] for j in range(len(power))]
        writer.writerow(row)

# you can plot directly here or use the .csv files in the coeff_plots.py script
# color blind friendly plotting color palette
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
plt.savefig('tests.pdf')
print('boop')
#plt.show()