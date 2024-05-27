import sys
import os
# Get the current directory of the script
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
hydro_dir = os.path.join(parent_dir, 'hydro')
swan_dir = os.path.join(parent_dir, 'swan')
sys.path.append(hydro_dir)
sys.path.append(swan_dir)

import sheep
import alpaca
import numpy as np

xgrid = 3000                                # size of grid in x-direction
ygrid = 5000                                # size of grid in y-direction
mxc = 300                                   # number of grid points in x (-1)
myc = 300                                   # number of grid points in y (-1)

x = [1490, 1590, 1690, 1420, 1520, 1620]    # x-position of bodies #, 1790,1890,1990] 
ya = 4500                                   # y-position of row 1
yb = 4400                                   # y-position of row 2
H = 0.8                                     # avg significant wave height [m] (i had 1.3832  for some reason?)
T = 5                                       # avg wave period [s] from buoy 44097
w = np.array([2*np.pi/T])   # wave frequency

breakwtr=False
point_absorber=False
oscillating_surge=False
attenuator=True
farm=False
controls=False

Kt_H, Kr_H, w_vals = sheep.wec_run(w,breakwtr,point_absorber,oscillating_surge,attenuator,farm,controls)
if farm == True:
    i = 1
else:
    i = 0

KR = [abs(Kr_H[0][0]), abs(Kr_H[i][0]), abs(Kr_H[2*i][0]), abs(Kr_H[0][0]), abs(Kr_H[i][0]), abs(Kr_H[2*i][0])]
KT = [abs(Kt_H[0][0]), abs(Kt_H[i][0]), abs(Kt_H[2*i][0]), abs(Kt_H[0][0]), abs(Kt_H[i][0]), abs(Kt_H[2*i][0])]
print(KR)
print(KT)

waveHeight = alpaca.swanrun(KR, KT, x, ya, yb, H, T, xgrid, ygrid, mxc, myc)