import numpy as np
import matplotlib.pyplot as plt
import PA
import OSWEC
import symsol
import antisol
import warnings
warnings.filterwarnings("ignore", category=np.ComplexWarning)

# step 1: import hydro coeffs, RAO, and radiation potential
w = 0.9             # wave frequency [rad/s]
dataset, rad_result, rad_prob, diff_prob, diff_result, body, pitch_RAO = OSWEC.bodysolver(w)
e = OSWEC.exc_force(diff_prob, rad_result, body, w)
A, B, K, M = OSWEC.hydro(dataset, body)
phi_rad, lam, grid = OSWEC.phi_rad(w, rad_result, pitch_RAO)

############# PROBLEM AREA AGAIN. WHERE TO PICK YOUR POINT ##################################
t, r, R, T = antisol.coefficients(phi_rad, B, w, M, A, K, e, lam)

plt.plot(np.linspace(0,4*lam,2*lam-1),t,color='black')
plt.plot(np.linspace(0,4*lam,2*lam-1),r,color='red')
plt.plot(np.linspace(0,4*lam,2*lam-1),T,color='black', linestyle=':')
plt.plot(np.linspace(0,4*lam,2*lam-1),R,color='red', linestyle=':')
plt.show()

#print('A_down',abs(A_down))
#print('hydrodynamic efficiency', 1 - abs(R)**2 - abs(T)**2)
#print('sum squared oscil', abs(R + T)**2)