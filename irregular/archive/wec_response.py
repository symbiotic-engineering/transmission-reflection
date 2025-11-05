# CALCULATING THE RAO FOR SPECTRAL BODY RESPONSE
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

breakwtr, point_absorber, oscillating_surge, attenuator = False, True, False, False
farm, controls, staggered, reactive = False, False, False, False
w = np.linspace(0.2,4,30)
S_j = jonswap.get_spectra(w)
RAO = []
omega = []
freq = []

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
    f = w / (np.pi * 2)
    RAO.append(RAO_vals[0][0])
    omega.append(w)
    freq.append(f)

irr_response = []

for i in range(np.size(S_j)):
    S_R = (RAO[i])**2 * S_j[i]
    irr_response.append(S_R)

# expected amplitude
xi = np.sqrt(simps(irr_response, omega))
print('expected body amplitude', np.abs(xi))

#plt.plot(omega,RAO)
plt.figure()
plt.plot(freq,S_j,label='JONSWAP')
plt.plot(freq, np.abs(irr_response), label='Body Response')
plt.legend()
plt.xlabel('f [Hz]')
plt.ylabel('Magnitude [m^2/Hz]')
plt.savefig('irr_resp.pdf')