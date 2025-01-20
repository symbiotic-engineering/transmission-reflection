# DISCLAIMER: this script needs a lot of work to be automated.
# currently, it's purpose is to give Maha preliminary irregular
# wave results to present at Georgia Tech
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
hydro_dir = os.path.join(parent_dir, 'hydro')
sys.path.append(hydro_dir)
import body
import solve
import jonswap

breakwtr, point_absorber, oscillating_surge, attenuator = False, False, True, False
farm, controls, staggered, reactive = False, True, False, True
w = np.array([0.7,0.8,0.9,1.1,1.25,1.3]) # 0.7,0.8,0.9,1.0, behavior a little weird at high freq, removing for now
omega = []
freq = []
H_ref = []
H_inc_up = []
H_trans = []
H_inc_down = []

# copy and pasted from "run_coeffs.py" in order to save time by not computing the free
# surface mesh. will be condensed in future iterations
N = 3
if farm == False:
    index = 1
else:
    index = 3
B = 0                                           # wave direction [rad]
depth = 500                                     # keep deep water assumption for EB
if staggered:
    xtrans = np.array([50,50])                  # x translation of bodies if staggered farm
    x_center = -25
else:
    xtrans = np.array([0,0])
    x_center = 0
ytrans = np.array([50,-50])                     # y translation of bodies if farm
for w in w:
    if w < 1.0:
        res = 2.0
    else:
        res = 3.5
    if breakwtr:
        array, rel_dim, char_dim = body.breakwater(xtrans,ytrans,farm,x_center)
        rad = False               # rad only false for breakwater case
    else:
        rad = True
    if point_absorber:
        array, rel_dim, char_dim, budal_limit = body.PA(xtrans,ytrans,farm,w,x_center)
    if oscillating_surge:
        array, rel_dim, char_dim, budal_limit = body.OSWEC(xtrans,ytrans,farm,w,x_center)
    if attenuator:
        array, rel_dim, char_dim, budal_limit = body.attenuator(xtrans,ytrans,farm,w,x_center)
    diff_result,rad_result,RAO_vals,lam,CWR = solve.hydro(array,B,depth,w,char_dim,farm,controls,point_absorber,reactive)
    total,incoming_fse,x1,x2,nx,y1,y2,ny = solve.elevation(res,lam,diff_result,rad_result,RAO_vals,farm,rad,controls,N,attenuator,rel_dim)

    import warnings
    warnings.filterwarnings("ignore", category=np.ComplexWarning)
    
    ##################################################################
    # Extract relevant columns
    mid_y = int(ny / 2)
    mid_x = int(nx / 2)
    convx = int(nx/(abs(x1)+x2))      # to convert meters to grid points
    convy = int(ny/(abs(y1)+y2))

    g = 9.81                # gravitational constant (m/s^2)
    k = w**2/g              # wave number infinite depth (rad^2/m)
    lam = int(2*np.pi/k)    # wavelength infinite depth (m)

    zinc_up = incoming_fse[mid_y, (mid_x - int(rel_dim*convx) + int(x_center*convx)) - int(lam*convx):mid_x - int(rel_dim*convx) + int(x_center*convx)]      # incident wave height upstream
    zinc_down = incoming_fse[mid_y, mid_x + int(rel_dim*convx) + int(x_center*convx):(mid_x + int(rel_dim*convx) + int(x_center*convx)) + int(lam*convx)]    # incident wave height downstream
    z_up = total[mid_y, (mid_x - int(rel_dim*convx) + int(x_center*convx)) - int(lam*convx):mid_x - int(rel_dim*convx) + int(x_center*convx)]                # total wave height upstream
    z_down = total[mid_y, mid_x + int(rel_dim*convx) + int(x_center*convx):(mid_x + int(rel_dim*convx) + int(x_center*convx)) + int(lam*convx)]              # transmitted wave height

    ''' this needs to somehow be a probability distribution function so i can get 
    the probability of a sea state occurring. since the JONSWAP spectra has
    units attached [m^2/Hz] idk how that works '''

    # S_j = jonswap.get_spectra(w)
    # S_H_ref = (z_up - zinc_up) * np.sqrt(S_j)
    # S_H_up = (zinc_up) * np.sqrt(S_j)
    # S_H_trans = (z_down) * np.sqrt(S_j)
    # S_H_down = (zinc_down) * np.sqrt(S_j)

    avg_H_zup = np.array([np.mean(abs(z_up))]) # - S_H_up))])
    avg_H_zincup = np.array([np.mean(abs(zinc_up))])
    avg_H_zdown = np.array([np.mean(abs(z_down))])
    avg_H_zincdown = np.array([np.mean(abs(zinc_down))])
    print('avg ref',avg_H_zup)
    print('avg inc up',avg_H_zincup)
    print('avg trans',avg_H_zdown)
    print('avg inc down',avg_H_zincdown)
    
    H_ref.append(avg_H_zup[0])
    H_inc_up.append(avg_H_zincup[0])
    H_trans.append(avg_H_zdown[0])
    H_inc_down.append(avg_H_zincdown[0])
    
    f = w / (np.pi * 2)
    omega.append(w)
    freq.append(f)

# expected wave heights
exp_H_ref = np.sqrt(simps(H_ref, freq))
exp_H_up = np.sqrt(simps(H_inc_up, freq))
exp_H_trans = np.sqrt(simps(H_trans, freq))
exp_H_down = np.sqrt(simps(H_inc_down, freq))
print('expected H_ref', exp_H_ref)
print('expected H_up',exp_H_up)
print('expected H_trans',exp_H_trans)
print('expected h_down',exp_H_down)

# expected coefficients
exp_kt = exp_H_trans/exp_H_down
exp_kr = (exp_H_ref/exp_H_up) - 1
print('expected transmission coefficient',exp_kt)
print('expected reflection coefficient',exp_kr)

# irregular coeffs (across the spectrum for plotting purposes)
irr_kt = []
irr_kr = []
for i in range(np.size(H_trans)):
    irr_KT = H_trans[i]/H_inc_down[i]
    irr_KR = H_ref[i]/H_inc_up[i]

    irr_kt.append(irr_KT)
    irr_kr.append(irr_KR)

colors = ['#377eb8', '#4daf4a']

plt.figure(figsize=(8, 6))
plt.plot(omega, np.abs(H_ref/H_inc_up) - 1, label='$K_r$', color=colors[0], linewidth=2)
plt.plot(omega, np.abs(H_trans/H_inc_down), label='$K_t$', color=colors[1], linewidth=2)
plt.legend(fontsize=18)
plt.xlabel('$\\omega$ [rad/s]', fontsize=20)
plt.ylabel('Wave Height', fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('irr_coeffs.pdf')