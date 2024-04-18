import sys
import os
import SWELLPA
import numpy as np
import matplotlib.pyplot as plt

# modeling WECs 1, 2, and 3 in SWELL dataset in RSS5 sea state

T = 1.5                 # period in seconds
w = np.array([2*np.pi/T])   # wave frequency
xtrans = np.array([-0.39,-0.78])
ytrans = np.array([0,0])

for w in w:
    if w < 0.8:
        res = 2
    else: 
        res = 4                                   # resolution factor of grid wrt lambda
    kd, total, incoming_fse, lam = SWELLPA.lpf(w,res,xtrans,ytrans)