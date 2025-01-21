import sys
import os
import SWELLPA
import numpy as np

# modeling WECs 1, 2, and 3 in SWELL dataset in RSS5 sea state
T = 1.5                                        # period [s]
w = np.array([2*np.pi/T])                      # wave frequency [rad/s]
xtrans = np.array([-0.39,-0.78])               # x-locs of WECs 1 and 3
ytrans = np.array([0,0])                       # y-locs of WECs 1 and 3
depth = 50                                     # water depth [m] 0.9

for w in w:
    total, incoming_fse, lam, elevation_at_gauges = SWELLPA.lpf(w,xtrans,ytrans,depth,farm=False)