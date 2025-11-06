import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from scipy.integrate import simps
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
hydro_dir = os.path.join(parent_dir, 'hydro')
sys.path.append(hydro_dir)
import body
import solve
import piersonmos
import pandas as pd
import csv

w = np.array([0.52359878,0.57119866,0.62831853,0.6981317,0.78539816,0.8975979,1.04719755,1.25663706])  # wave frequencies corresponding to SouthFork spectrum
def get_coefficients(csv_file):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_file)
    exp_Kt = [[],[],[],[]]
    exp_Kr = [[],[],[],[]]
    KT = []
    KR = []
    S_KR = []
    S_KT = []

    S_pm = piersonmos.get_spectra(w)

    for j in range(4):
        S_Kr = [[],[],[],[],[],[],[],[]]
        S_Kt = [[],[],[],[],[],[],[],[]]
        Kt = [[],[],[],[],[],[],[],[]]
        Kr = [[],[],[],[],[],[],[],[]]
        for i in range(8):
            Kt[i] = float(df.iloc[7-i, j+1])
            Kr[i] = float(df.iloc[7-i, j+5])
            S_Kr[i] = S_pm[i] * Kr[i]**2
            S_Kt[i] = S_pm[i] * Kt[i]**2

        # expected coefficients
        exp_Kr[j] = sum(cumulative_trapezoid(S_Kr, w))
        exp_Kt[j] = sum(cumulative_trapezoid(S_Kt, w))

        KT.append(Kt)
        KR.append(Kr)
        S_KR.append(S_Kr)
        S_KT.append(S_Kt)

    print('expected Kr', exp_Kr)
    print('expected Kt',exp_Kt)

    return exp_Kt, exp_Kr, KT, KR, S_KR, S_KT

csv_file_name = 'PA_spectra_damp.csv'
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
hydro_data_dir = os.path.join(base_dir, 'hydro', 'data')
csv_file = os.path.join(hydro_data_dir, csv_file_name)

data_folder = os.path.join(os.path.dirname(__file__), 'data')
file_name = csv_file_name
file_path = os.path.join(data_folder, file_name)

exp_Kt, exp_Kr, KT, KR, S_KR, S_KT = get_coefficients(csv_file)

wS = np.linspace(0.45,1.3,40)
S_pm = piersonmos.get_spectra(wS)


colors = ['#377eb8', '#4daf4a']

plt.figure(figsize=(8, 6))
plt.plot(w, S_KT[0], marker = 'p', label='$K_t$ WEC1', color=colors[0], linewidth=2)
plt.plot(w, S_KR[0], marker = 'p', label='$K_r$ WEC2', color=colors[1], linewidth=2)
plt.plot(w, S_KT[1], marker = 'o', label='WEC2', color=colors[0], linewidth=2)
plt.plot(w, S_KR[1], marker = 'o', color=colors[1], linewidth=2)
plt.plot(w, S_KT[2], marker = 's',label='WEC3', color=colors[0], linewidth=2)
plt.plot(w, S_KR[2], marker = 's', color=colors[1], linewidth=2)
plt.plot(w, S_KT[3], marker = '>', label='WEC4', color=colors[0], linewidth=2)
plt.plot(w, S_KR[3], marker = '>', color=colors[1], linewidth=2)
plt.plot(wS,S_pm,label='PM Spectrum',color='black',linewidth=2)
plt.legend(fontsize=18)
plt.xlabel('$\omega$ [rad/s]', fontsize=20)
plt.ylabel('Spectrum [$m^2s/rad$]', fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('irr_coeffs.pdf')

# Save the data to the corresponding .csv file
with open(file_path, mode='w', newline='') as file:
    writer = csv.writer(file)
    header = ['Device Number'] + [f'Kt_{i+1}' for i in range(len(exp_Kt))] + [f'Kr_{i+1}' for i in range(len(exp_Kr))]
    writer.writerow(header)
    row = exp_Kt + exp_Kr
    writer.writerow(row)