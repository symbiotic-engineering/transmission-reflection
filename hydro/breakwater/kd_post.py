import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

w = 1.047
g = 9.81
k = w**2 / g
lambda_wave = 2 * np.pi / k
int_lam = int(lambda_wave)
xR = np.linspace(-1*int(int_lam),0,int(int_lam // 2))
#print(xR)
xD = np.linspace(0,1*int(int_lam),int(int_lam // 2))
xt = np.linspace(-1*int(int_lam),1*int(int_lam),int(int_lam))
y0 = np.zeros_like(xt)

#########################################################################################################
################ RANDOM TEST FOR KD #######################
# Load datasets
inc_el = "incoming_fse_data_bw.csv"
disturb = "kd_1047.csv"
inc_el = pd.read_csv(inc_el, header=None)
kd1 = pd.read_csv(disturb, header=None)

# Extract relevant columns
zinc_up = inc_el.applymap(lambda x: complex(x.replace(' ', ''))).iloc[:int(int_lam // 2),int(int_lam // 2)-1]
kd_up = kd1.applymap(lambda x: complex(x.replace(' ', ''))).iloc[:int(int_lam // 2),int(int_lam // 2)-1]

def find_peaks(data, lambda_wave, xR):
    peaks = []
    for i in range(len(data)):
        if i < 0 or i >= len(data):
            continue
        if i == 0 or i == len(data) - 1:
            continue
        if xR[i] < -1*int(lambda_wave)//2 and data[i] > data[i - 1] and data[i] > data[i + 1]:
            peaks.append(data[i])
    return peaks
def avg_peak_height(data, lambda_wave, xR):
    peak_values = find_peaks(data, lambda_wave, xR)
    average_peak_height = np.mean(peak_values)
    return average_peak_height


avg_kd_up = avg_peak_height(kd_up, lambda_wave, xR)
avg_zinc_up = avg_peak_height(zinc_up, lambda_wave, xR)

reflection_coeff = abs(avg_kd_up)
print('ref coef',reflection_coeff)

######################################## FOR DOWNSTREAM ###########################################
zinc_down = inc_el.applymap(lambda x: complex(x.replace(' ', ''))).iloc[int(int_lam // 2):,int(int_lam // 2)-1]
kd_down = kd1.applymap(lambda x: complex(x.replace(' ', ''))).iloc[int(int_lam // 2):,int(int_lam // 2)-1]

def peaks_down(series, lambda_wave):
    peaks = []
    for i in range(1, len(series) - 1):  
        if np.angle(series[i] - series[i - 1]) > 0 and np.angle(series[i] - series[i + 1]) > 0:
            if xD[i] < lambda_wave:  # Limiting to x < lambda_wave
                peaks.append(series[i].real)  
    return peaks

def avg_peak_down(data, lambda_wave):
    peak_values = peaks_down(data, lambda_wave)
    avg_peak_values = np.mean(peak_values)
    return avg_peak_values

avg_distance_kd_down = avg_peak_down(kd_down.values, lambda_wave)
avg_distance_zinc_down = avg_peak_down(zinc_down.values, lambda_wave)

trans_coeff = abs(avg_distance_kd_down)
print('trans coef',trans_coeff)

EB1 = trans_coeff**2 + (reflection_coeff - 1)**2       # energy balance check
EB2 = (trans_coeff + (reflection_coeff - 1))**2       # energy balance check
print('energy balance, sum of squares',EB1)
print('energy balance, sum squared',EB2)

lam_len_neg = np.linspace(-int_lam*1,0)
lam_len_pos = np.linspace(0,int_lam*1)

avg_up = reflection_coeff * np.ones_like(lam_len_neg)
avg_down = trans_coeff * np.ones_like(lam_len_pos)

plt.plot(xR, kd_up, color='orange', label='disturbance coeff')
plt.plot(xD, kd_down, color='orange')
plt.plot(lam_len_neg,avg_up,color='black',linestyle=':',label='Peak Avg')
plt.plot(lam_len_pos,avg_down,color='black',linestyle=':')
plt.plot(xR, zinc_up, color='blue', linestyle=':', label='incident fse')
plt.plot(xD, zinc_down, color='blue', linestyle=':')
# plt.plot(xt, y0, color='green', linestyle='--', label='Still Water Line')
plt.xlabel('x [m]')
plt.ylabel('Elevation')
plt.legend()
plt.savefig(f'kd/1047_kdcomp.pdf')
plt.show()

#################################### PEAKS AND TROUGHS ###################################
# # finds peaks and troughs
# def find_peaks_and_troughs(data, range_indices):
#     peaks = []
#     troughs = []
#     for i in range_indices:
#         if i < 0 or i >= len(data):
#             continue
#         if i == 0 or i == len(data) - 1:
#             continue
#         if data[i] > data[i - 1] and data[i] > data[i + 1]:
#             peaks.append(data[i])
#         elif data[i] < data[i - 1] and data[i] < data[i + 1]:
#             troughs.append(data[i])
#     return peaks, troughs

# def avg_peak_and_trough_values(data, lambda_wave):
#     # Find indices within ±2λ range
#     center_index = len(data)
#     range_indices = range(center_index - int(lambda_wave), center_index)
#     # print(range_indices)

#     # Find peaks and troughs within the specified range
#     peaks_values, troughs_values = find_peaks_and_troughs(data, range_indices)

#     # Calculate the average of peak and trough values
#     avg_peak_trough_values = np.mean(peaks_values + troughs_values)
#     return avg_peak_trough_values

# # Calculate average wave heights for each dataset
# avg_kd_up = avg_peak_and_trough_values(kd_up, lambda_wave)
# avg_zinc_up = avg_peak_and_trough_values(zinc_up, lambda_wave)

# reflection_coeff = abs(avg_kd_up)
# print('ref coef',reflection_coeff)

# # Define a function to find peaks and troughs within a given range
# def find_peaks_and_troughs(series, lambda_wave):
#     peaks = []
#     troughs = []
#     for i in range(1, len(series) - 1):  # Start from index 1 and end at len(series) - 2
#         if np.angle(series[i] - series[i - 1]) > 0 and np.angle(series[i] - series[i + 1]) > 0:
#             if i < lambda_wave:  # Limiting to x < 2 * lambda
#                 peaks.append(series[i].real)  # Append only the real part of the peak
#         elif np.angle(series[i] - series[i - 1]) < 0 and np.angle(series[i] - series[i + 1]) < 0:
#             troughs.append(series[i].real)  # Append only the real part of the trough
#     return peaks, troughs

# def avg_peak_and_trough_values(data, lambda_wave):
#     # Find peaks and troughs within the specified range
#     peaks_values, troughs_values = find_peaks_and_troughs(data, lambda_wave)

#     # Calculate the average of peak and trough values
#     avg_peak_trough_values = np.mean(peaks_values + troughs_values)
#     return avg_peak_trough_values

# # Calculate average of peak and trough values for each dataset
# avg_distance_kd_down = avg_peak_and_trough_values(kd_down.values, lambda_wave)
# avg_distance_zinc_down = avg_peak_and_trough_values(zinc_down.values, lambda_wave)
