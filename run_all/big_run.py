''' TO RUN THIS SCRIPT (read the README too): 
    1. You need to define your xgrid and ygrid size and resolution.
    2. You need to check the xtrans and ytrans WEC positions in the
    sheep.wec_run() function. This will determine if you will run a 
    staggered or regular array.
    3. You need to define the positions of your bodies in the SWAN grid
    according to the sheep.wec_run() function.
    4. You need to define the significant wave height and peak wave
    period of your sea state. ATTN: you will also need to adjust the
    run_swan.py file for wind speed and direction.
    5. Decide which body you are running by setting it equal to True.
    Set all others to False.
    6. If you want to include array interactions, set farm equal to
    True. If you want to include controls, set controls equal to True.
    ATTN: whether you use pure damping or reactive controls is determined
    by the PTO.py script in the hydro folder.'''

import sys
import os
import run_coeffs
import post_process
import run_swan
import numpy as np

# Get the current directory of the script
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
hydro_dir = os.path.join(parent_dir, 'hydro')
swan_dir = os.path.join(parent_dir, 'swan')
sys.path.append(hydro_dir)
sys.path.append(swan_dir)

xgrid = 3000                                # size of grid in x-direction
ygrid = 5000                                # size of grid in y-direction
mxc = 300                                   # number of grid points in x (-1) (10 m res in x)
myc = 500                                   # number of grid points in y (-1) (10 m res in y)

# six-body reg: 
# x = [1290,1390,1490,1590,1690,1790] # x-position of bodies
# ya=4500         # y-position of row 1
# yb=4500         # y-position of row 2

# three-body stag: x = [1490,1590,1690]ya=4500 yb=4450 # change i==1 in run_swan

# three-body reg: x=[1490,1590,1690] ya=4500 yb=4400

# six body stag: 
x=[1490, 1590, 1690, 1420, 1520, 1620]
ya=4500 
yb=4400
                                  
H = 0.8                                     # avg significant wave height [m]
T = 5                                       # avg wave period [s] from buoy 44097
w = np.array([2*np.pi/T])   # wave frequency

breakwtr=False
point_absorber=False
oscillating_surge=True
attenuator=False

farm=True
controls=True

Kt_H, Kr_H, w_vals = run_coeffs.wec_run(w,breakwtr,point_absorber,oscillating_surge,attenuator,farm,controls)
if farm == True:
    i = 1
else:
    i = 0

KR = [Kr_H[0][0], Kr_H[i][0], Kr_H[2*i][0], Kr_H[0][0], Kr_H[i][0], Kr_H[2*i][0]]
KT = [Kt_H[0][0], Kt_H[i][0], Kt_H[2*i][0], Kt_H[0][0], Kt_H[i][0], Kt_H[2*i][0]]
print(KR)
print(KT)

# diameter of the bodies
if point_absorber:
    d = 21
if oscillating_surge:
    d = 18
if breakwtr:
    d = 20
if attenuator:
    d = 4

sfgrid_dat, sfgrid_tbl = run_swan.generate_swan_input(KR, KT, d, x, ya, yb, H, T, xgrid, ygrid, mxc, myc,attenuator)
waveHeight = post_process.postpro(sfgrid_dat,xgrid,ygrid,mxc,myc)