import numpy as np
import matplotlib.pyplot as plt
import PA
import OSWEC
import attenuator
import symsol
import antisol
import warnings
warnings.filterwarnings("ignore", category=np.ComplexWarning)

# step 1: import hydro coeffs, RAO, and radiation potential
w = 1.047             # wave frequency [rad/s]
dataset, rad_result, rad_prob, diff_prob, diff_result, body = PA.bodysolver(w)
e = PA.exc_force(diff_prob, rad_result, body, w)
A, B, K, M = PA.hydro(dataset, body)
phi_rad, lam, grid = PA.phi_rad(w, rad_result, diff_result)

############# PROBLEM AREA AGAIN. WHERE TO PICK YOUR POINT ##################################
t, r, R, T = symsol.coefficients(phi_rad, B, w, M, A, K, e, lam)

plt.plot(np.linspace(0,4*lam,2*lam-1),t,color='black')
plt.plot(np.linspace(0,-4*lam,2*lam-1),r,color='red')
plt.plot(np.linspace(0,4*lam,2*lam-1),T,color='black', linestyle=':')
plt.plot(np.linspace(0,-4*lam,2*lam-1),R,color='red', linestyle=':')
plt.show()

avg_R, avg_T, avg_r, avg_t = symsol.avg_coeff(t,r,R,T,lam)
print('R',avg_R)
print('T',avg_T)
print('r',avg_r)
print('t',avg_t)

#print('A_down',abs(A_down))
#print('hydrodynamic efficiency', 1 - abs(R)**2 - abs(T)**2)
#print('sum squared oscil', abs(R + T)**2)