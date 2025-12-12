''' this determination of sliding friciton in the PA linear sleeve bearings is adapted
from the Nature paper by Vitorino 2017 titled:
'Effect of Sliding Friction in Harmonic Oscillators'
'''
'''
FURTHER WORK WITH THIS MODEL DEEMED IT WRONG. THE FRICTION DAMPING COEFF
MUST BE LARGER WITH FASTER FREQUENCIES, AND THAT IS NOT THE OBSERVED
RELATIONSHIP FROM THIS MODEL
'''
def slidefriction(B,M,K,A,w,surge_force):
    import numpy as np

    F_f = surge_force           # initial friction force estimation
    n = np.array([1,3]) # fourier series term number
    x_0 = 1             # excitation amplitude --> set to 1 for capy unit amplitude convention

    gam_n = []
    for n in n:
        gam_loop = ((4*F_f)/np.pi) * ( ((M + A) * n**2 * w**2 - K) / (B**2 * n**2 * w**2 + (K - (M + A) * n**2 * w**2)**2) )
        gam_n.append(gam_loop)
    
    gam = np.sum(gam_n)

    Z = np.sqrt( (K - (M + A) * w**2)**2 + B**2 * w**2 )

    A_n = ( -4 * F_f * B * w + np.sqrt(np.pi**2 * K**2 * x_0**2 * Z**2 - (4 * F_f * (K - (M + A) * w**2) - np.pi * gam * Z**2)**2) ) / (np.pi * Z**2)

    R = np.sqrt(A_n**2 + gam**2)

    B_slide = (4 * F_f) / (np.pi * w * R)

    return B_slide

