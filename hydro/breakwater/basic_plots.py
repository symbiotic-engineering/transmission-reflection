import matplotlib.pyplot as plt
import numpy as np

# plotting collected Kt and Kr values
# i am aware this is not efficient but i need to find this trend
Kt_wh = np.array([1.0089, 0.5913, 0.3752, 0.1872, 0.1361])   # 0.1712, 0.1361])        # from wave height
Kr_wh = np.array([0.0747, 0.4496, 0.6352, 0.7605, 0.8817])    # 1.2797, 0.8817])        # from wave height

Kt_kd = np.array([0.9465, 0.2835, 0.2013, 0.2467, 0.1247])  # 0.1704, 0.1247])        # from disturbance
Kr_kd = np.array([0.2062, 0.5067, 0.6490, 0.7524, 0.8475])  #1.3222, 0.8475])       # from disturbance

w = np.array([0.5, 0.65, 0.75, 0.85, 1.047])    # 0.95, 1.047])                        # frequencies

plt.plot(w,Kt_wh,color='blue',marker='o',label='Kt from H')
plt.plot(w,Kr_wh,color='red',marker='o',label='Kr from H')
plt.plot(w,Kt_kd,color='blue',marker='o',linestyle=':',label='Kt from Kd')
plt.plot(w,Kr_kd,color='red',marker='o',linestyle=':',label='Kr from Kd')
plt.legend()
plt.xlabel('w [rad/s]')
plt.ylabel('Coefficient Value [m]')
plt.savefig(f'plots/trans_ref_coeffs_yay.pdf')
plt.show()
