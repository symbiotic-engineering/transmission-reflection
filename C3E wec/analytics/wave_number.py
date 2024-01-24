# finding wave number and wavelength for infinite and finite depth (kH --> 0)
import numpy as np

w = np.array([0.5,0.75,0.85,0.9,1.047,1.2,1.35,1.5,1.75,2.0])               #np.array([0.5, 1, 1.047, 1.5, 2, 2.5, 3])
H = 40                                #40

def wave_num(w,H):
    # for infinite depth
    g = 9.81            # gravitational constant (m/s^2)
    k_inf = w**2/g      # wave number infinite depth (rad^2/m)
    lam_inf = 2*np.pi/k_inf    # wavelength infinite depth (m)

    # for finite depth
    k_f = w/np.sqrt(g*H)    # wave number in finite depth
    lam_f = 2*np.pi/k_f     # wave length in finite depth
    return[k_inf,lam_inf,k_f,lam_f]



[kinf, laminf, kf, lamf] = wave_num(w,H)

print('infinite depth k', kinf)
print('infinite depth lamda', laminf)
print('finite depth k', kf)
print('finite depth lamda',lamf)