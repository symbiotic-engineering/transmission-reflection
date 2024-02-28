import matplotlib.pyplot as plt
import numpy as np

# plotting collected Kt and Kr values
# i am aware this is not efficient but i need to find this trend
Kt_wh = np.array([0.897053, 0.826340, 0.6833328, 0.689504, 0.4812433, 0.5166075])      # from wave height
Kr_wh = np.array([0.001443, 0.019126, 0.0674371, 0.2523389, 0.7058256, 0.70128805])      # from wave height
#EB1 []
#EB2 []

Kt_kd = np.array([0.992879, 0.909906, 0.786848, 0.128409, 0.0138686, 0.3114734])  # from disturbance
Kr_kd = np.array([0.024969, 0.081997, 0.172210, 0.345168, 0.5896767, 0.7011844])  # from disturbance
#EB1 []
#EB2 []

w = np.array([0.5, 0.65, 0.75, 0.85, 0.95, 1.047])      # frequencies

plt.plot(w,Kt_wh,color='blue',marker='o',label='Kt from H')
plt.plot(w,Kr_wh,color='red',marker='o',label='Kr from H')
plt.plot(w,Kt_kd,color='blue',marker='o',linestyle=':',label='Kt from Kd')
plt.plot(w,Kr_kd,color='red',marker='o',linestyle=':',label='Kr from Kd')
plt.legend()
plt.xlabel('w [rad/s]')
plt.ylabel('Coefficient Value [m]')
plt.savefig('trans_ref_coeffs_yay.pdf')
plt.show()
