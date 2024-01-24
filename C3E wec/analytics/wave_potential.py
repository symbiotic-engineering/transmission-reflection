# finding the wave potential in finite and infinite depth
# obtain k_inf and k_f from wave_num
# w: omega [rad/s]
# t: time [s]
# H: ocean depth
# A: incoming wave amplitude

import numpy as np
import matplotlib.pyplot as plt

#w = np.array([0.5, 1, 1.047, 1.5, 2, 2.5, 3])
w = 1.047
H = 40

def wave_num(w,H):
    # for infinite depth
    g = 9.81            # gravitational constant (m/s^2)
    k_inf = w**2/g      # wave number infinite depth (rad^2/m)
    lam_inf = 2*np.pi/k_inf    # wavelength infinite depth (m)

    # for finite depth
    k_f = w/np.sqrt(g*H)    # wave number in finite depth
    lam_f = 2*np.pi/k_f     # wave length in finite depth
    return[k_inf,lam_inf,k_f,lam_f]

[k_inf,lam_inf,k_f,lam_f] = wave_num(w,H)

t = np.linspace(0,20,num=40)
A = 1
x = 1
z = -20

def potential(k_inf,k_f,w,t,x,z,H,A):
    # for infinite depth
    g = 9.81                # gravitational constant [m/s^2]
    phi_inf_complex = ((g*A*1j)/w)*np.exp(k_inf*z - k_inf*x*1j + w*t*1j)
    phi_inf_real = phi_inf_complex.real

    # for finite depth
    phi_f_complex = ((g*A*1j)/w)*((np.cosh(k_f*(z+H)))/(np.cosh(k_f*H)))*np.exp(-k_f*x*1j + w*t*1j)
    phi_f_real = phi_f_complex.real
    return(phi_inf_real,phi_f_real)

[p_inf,p_f] = potential(k_inf,k_f,w,t,x,z,H,A)

# plot potential
plt.plot(t,p_inf,'r', label = 'Infinite Depth')
plt.plot(t,p_f,'k', label = 'Finite Depth')
plt.xlabel('t [s]',fontsize=17)
plt.ylabel('Potential',fontsize=17)
plt.legend()
plt.show()
