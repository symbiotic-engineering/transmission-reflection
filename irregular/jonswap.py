# IRREGULAR WAVE MODELING FOR NEAR FIELD
import numpy as np
import matplotlib.pyplot as plt

def get_spectra(w):
    g = 9.81                                    # gravitational constant [m/s^2]
    U_10 = 4.92                                 # wind speed [m/s]
    F = 104956.99                               # dist where wind has const velocity [m/s] calculated from peak freq eq
    alpha = 0.076 * (U_10**2 / (F * g))**0.22   # factor in spectral eq
    w_p = (2*np.pi)/5                           # peak frequency [rad/s]
    gamma = 1.54                                # enhancement factor
    S_j = []                                    # JONSWAP wave spectral density [m^2/Hz]

    for omega in w:
        if omega <= w_p:
            sigma = 0.07
        else:
            sigma = 0.09
        r = np.exp(-(omega - w_p)**2 / (2 * sigma**2 * w_p**2))

        S = ((alpha*g**2)/omega**5) * np.exp((-5/4) * (w_p/omega)**4) * gamma**r 
        S_j.append(S)

    f = w / (np.pi * 2)
    plt.plot(f,S_j)
    plt.savefig('jonswap.pdf')

    return S_j