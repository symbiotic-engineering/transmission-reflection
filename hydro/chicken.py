# for comparing controlled and uncontrolled WECs
import numpy as np
import matplotlib.pyplot as plt
import sheep

file_path = '/mnt/c/Users/ov162/transmission-reflection/hydro/figures/'
file_name = 'atten_cont_comp.pdf'

breakwtr=False
point_absorber=False
oscillating_surge=False
attenuator=True
farm=False

w = np.array([0.5,0.65,0.75,0.85,0.95,1.047,1.2])   # wave frequency

Kt_H, Kr_H, w_vals = sheep.wec_run(w,breakwtr,point_absorber,oscillating_surge,attenuator,farm,controls=False)
Kt_cont, Kr_cont, w_vals = sheep.wec_run(w,breakwtr,point_absorber,oscillating_surge,attenuator,farm,controls=True)

cud_colors = ['#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7', '#000000', '#8B4513']
linestyles = ['-', '--', ':', '-.', '-', '--', ':', '-.']

for i, kt_h_values in enumerate(Kt_H):
    plt.plot(w_vals, kt_h_values, marker='o', label=f'$K_t$ body {i + 1}', 
             color=cud_colors[i % len(cud_colors)], linestyle=linestyles[i % len(linestyles)])
for i, kr_h_values in enumerate(Kr_H):
    plt.plot(w_vals, kr_h_values, marker='x', label=f'$K_r$ body {i + 1}', 
             color=cud_colors[i+1 % len(cud_colors)], linestyle=linestyles[i % len(linestyles)])
for i, kt_cont in enumerate(Kt_cont):
    plt.plot(w_vals, kt_cont, marker='*', label=f'Controlled $K_t$ body {i + 1}', 
             color=cud_colors[i+2 % len(cud_colors)], linestyle=linestyles[i % len(linestyles)])
for i, kr_cont in enumerate(Kr_cont):
    plt.plot(w_vals, kr_cont, marker='+', label=f'Controlled $K_r$ body {i + 1}', 
             color=cud_colors[i+3 % len(cud_colors)], linestyle=linestyles[i % len(linestyles)])


plt.legend()
plt.xlabel('$\omega$ [rad/s]')
plt.ylabel('Coefficient Value')
plt.savefig(f'{file_path}{file_name}')
plt.show()