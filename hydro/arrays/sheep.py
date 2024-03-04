# run all Kt and Kr calcs for any body here
import PA_array
import OSWEC_array
import wave_height
import kd_post
import numpy as np
import matplotlib.pyplot as plt

w = np.array([0.5,0.65,0.75,0.85,0.95,1.047])   # wave frequency
xtrans = np.array([25,25])
ytrans = np.array([50,-50])

Kr_H = []
Kt_H = []
Kr_K = []
Kt_K = []
w_vals = []
for w in w:
    if w < 0.8:
        res = 2
    else: 
        res = 2                                   # resolution factor of grid wrt lambda
    kd, total, incoming_fse, lam = OSWEC_array.lpf(w,res,xtrans,ytrans)
    ref_H, trans_H, EB1, EB2 = wave_height.wave_height(total, incoming_fse,lam, res)
    ref_K, trans_K, EB1, EB2 = kd_post.disturbance(kd, lam, res)
    Kr_H.append(ref_H)
    Kt_H.append(trans_H)
    Kr_K.append(ref_K)
    Kt_K.append(trans_K)
    w_vals.append(w)

plt.plot(w_vals,Kt_H,color='blue',marker='o',label='Kt from H')
plt.plot(w_vals,Kr_H,color='red',marker='o',label='Kr from H')
plt.plot(w_vals,Kt_K,color='blue',marker='o',linestyle=':',label='Kt from Kd')
plt.plot(w_vals,Kr_K,color='red',marker='o',linestyle=':',label='Kr from Kd')
plt.legend()
plt.xlabel('w [rad/s]')
plt.ylabel('Coefficient Value [m]')
plt.show()