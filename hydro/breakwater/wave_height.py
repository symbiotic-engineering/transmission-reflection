import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=np.ComplexWarning)
##################################################################

w = np.array([0.50])

def wave_num(w):
    # for infinite depth
    g = 9.81            # gravitational constant (m/s^2)
    k_inf = w**2/g      # wave number infinite depth (rad^2/m)
    lam_inf = 2*np.pi/k_inf    # wavelength infinite depth (m)
    return k_inf, lam_inf

k, lambda_wave = wave_num(w)
int_lam = int(lambda_wave[0])
xR = np.linspace(-2*int_lam,0,int(int_lam // 2))
xD = np.linspace(0,2*int_lam,int(int_lam // 2))
#xt = np.linspace(-2*int_lam,2*int_lam,int_lam)
#y0 = np.zeros_like(xt)

# Load your datasets
capy_elev = "dif_plus_inc_data_bw.csv"
inc_el = "incoming_fse_data_bw.csv"

el = pd.read_csv(capy_elev, header=None)
inc_el = pd.read_csv(inc_el, header=None)

# Extract relevant columns
z_up = el.applymap(lambda x: complex(x.replace(' ', ''))).iloc[:int(int_lam // 2),int(int_lam // 2)-1]
zinc_up = inc_el.applymap(lambda x: complex(x.replace(' ', ''))).iloc[:int(int_lam // 2),int(int_lam // 2)-1]
z_down = el.applymap(lambda x: complex(x.replace(' ', ''))).iloc[int(int_lam // 2):,int(int_lam // 2)-1]
z_down.reset_index(drop=True, inplace=True)
zinc_down = inc_el.applymap(lambda x: complex(x.replace(' ', ''))).iloc[int(int_lam // 2):,int(int_lam // 2)-1]
zinc_down.reset_index(drop=True, inplace=True)

# Calculate average wave heights for each dataset
pk_zup, _ = find_peaks(z_up.values, height = 0)         # peak
pk_zincup, _ = find_peaks(zinc_up.values,height = 0)    # peak
tr_zup, _ = find_peaks(-z_up.values, height = 0)        # trough
tr_zincup, _ = find_peaks(-zinc_up.values,height = 0)   # trough
pk_zdown, _ = find_peaks(z_down.values, height = 0)         # peak
pk_zincdown, _ = find_peaks(zinc_down.values,height = 0)    # peak
tr_zdown, _ = find_peaks(-z_down.values, height = 0)        # trough
tr_zincdown, _ = find_peaks(-zinc_down.values,height = 0)   # trough

def calculate_average_wave_height(data):
    peaks, _ = find_peaks(data, height=0)
    troughs, _ = find_peaks(-data, height=0)

    closest_diff = float('inf')
    closest_pair = None
    H_values = []

    for pk_index in peaks:
        for tr_index in troughs:
            diff = abs(data[pk_index] - data[tr_index])
            if diff < closest_diff:
                closest_diff = diff
                closest_pair = (pk_index, tr_index)

    if closest_pair:
        peak_index, trough_index = closest_pair
        H = data[peak_index] - data[trough_index]
        H_values.append(H)

    average_H = np.mean(H_values) if H_values else 0
    return average_H

avg_H_zup = calculate_average_wave_height(z_up.values)
avg_H_zincup = calculate_average_wave_height(zinc_up.values)
avg_H_zdown = calculate_average_wave_height(z_down.values)
avg_H_zincdown = calculate_average_wave_height(zinc_down.values)
ref_coeff = abs(avg_H_zup/avg_H_zincup) - 1
trans_coeff = abs(avg_H_zdown/avg_H_zincdown)
EB1 = trans_coeff**2 + ref_coeff**2     # energy balnce
EB2 = (trans_coeff + ref_coeff)**2      # energy balance
print('reflection coefficient',ref_coeff)
print('transmission coefficient',trans_coeff)
print('energy balance sum of squares',EB1)
print('energy balance square of sums',EB2)

plt.plot(xR, z_up, color ='black', label='Total')
plt.plot(xR[pk_zup], z_up[pk_zup], "o", color='red')
plt.plot(xR[tr_zup], z_up[tr_zup], "o", color='red')
plt.plot(xR, zinc_up, color ='blue', label='Incident')
plt.plot(xR[pk_zincup], zinc_up[pk_zincup], "o", color='red')
plt.plot(xR[tr_zincup], zinc_up[tr_zincup], "o", color='red')
plt.plot(xD, z_down, color ='black', label='Total')
plt.plot(xD[pk_zdown], z_down[pk_zdown], "o", color='red')
plt.plot(xD[tr_zdown], z_down[tr_zdown], "o", color='red')
plt.plot(xD, zinc_down, color ='blue')
plt.plot(xD[pk_zincdown], zinc_down[pk_zincdown], "o", color='red')
plt.plot(xD[tr_zincdown], zinc_down[tr_zincdown], "o", color='red')
plt.savefig(f'new_plt/050test_elev.pdf')
plt.show()
