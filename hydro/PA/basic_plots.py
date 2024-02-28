import matplotlib.pyplot as plt
import numpy as np

# plotting collected Kt and Kr values
# i am aware this is not efficient but i need to find this trend
Kt_wh = np.array([1.001698 , 1.005225, 1.010683, 1.01661 , 1.022065, 1.02327])      # from wave height
Kr_wh = np.array([-0.002482, -0.008944, -0.015262, -0.0187300, -0.02978,-0.03099])      # from wave height

Kt_kd = np.array([1.002368 , 1.005798, 1.002825, 0.993599, 0.989658, 0.966108])  # from disturbance
Kr_kd = np.array([0.002775 , 0.0091095, 0.016434, 0.025248, 0.028436, 0.020151])  # from disturbance

w = np.array([0.5, 0.65, 0.75, 0.85, 0.95, 1.047])                        # frequencies

plt.plot(w,Kt_wh,color='blue',marker='o',label='Kt from H')
plt.plot(w,Kr_wh,color='red',marker='o',label='Kr from H')
plt.plot(w,Kt_kd,color='blue',marker='o',linestyle=':',label='Kt from Kd')
plt.plot(w,Kr_kd,color='red',marker='o',linestyle=':',label='Kr from Kd')
plt.legend()
plt.xlabel('w [rad/s]')
plt.ylabel('Coefficient Value [m]')
plt.savefig(f'trans_ref_coeffs_yay.pdf')
plt.show()
