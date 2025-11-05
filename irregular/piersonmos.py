# IRREGULAR WAVE MODELING FOR NEAR FIELD
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

def get_spectra(w):

    g = 9.81                                     # acceleration due to gravity [m/s]
    Hs = 2.19                                    # significant wave height [m]
    Tp = 11                                      # dominant wave period [s]
    wp = (2*np.pi)/Tp                            # dominant wave frequency [rad/s]
    alpha = 8.1*10**-3                           # dimless coeff in equation
    U_10 = 4.92                                  # wind speed at 10 m above water surface [m/s]
    U_195 = 1.026*U_10                           # wind speed at 19.5 m above water surface [m/s]

    ### modified pierson-moskowitz (bretschneider)
    ### derived to only depend on peak period and significant wave height
    ### what SWAN will use when calling PM spectrum
    S_pm = (5/16)*(Hs**2)*(wp**4)*(w**-5)*np.exp((-5/4)*(wp/w)**4)
    #S_pm = (alpha*g**2 / w**5) * np.exp((-5/4)*(wp/w)**4)

    integral_spm = cumulative_trapezoid(S_pm, w)
    print('integral val',sum(integral_spm))

    return S_pm

# # Main script
# w = np.linspace(0.45,1.3,40)#np.array([0.48332195,0.57119866,0.6981317,0.8975979,1.25663706])  # Frequencies in rad/s

# # Adjust font sizes
# plt.rcParams.update({
#     'font.size': 24,  # Default font size
#     'axes.labelsize': 24,  # Axes label font size
#     'xtick.labelsize': 20,  # X-tick label font size
#     'ytick.labelsize': 20,  # Y-tick label font size
#     'legend.fontsize': 18,  # Legend font size
#     'figure.titlesize': 22,  # Figure title font size
# })

# # Use a colorblind-friendly color palette
# colors = ["#0072B2", "#D55E00", "#CC79A7", "#009E73"]  # Blue, Vermilion, Yellow, Green

# plt.figure(figsize=(12, 8))
# spectra = get_spectra(w)                                             # compute spectra
# plt.plot(w, spectra, color=colors[0], linewidth=4)

# # Plot settings
# plt.xlabel('$\omega$ [rad/s]')
# plt.ylabel('Spectral Density [$m^2/\omega$]')
# plt.grid()
# plt.savefig('pm_spectrum.pdf')

    ##### these parameters are for finding the JONSWAP spectrum. but we are in deep water, so we need
    ##### pierson-moskowitz spectrum instead
    # g = 9.81                                    # gravitational constant [m/s^2]
    # U_10 = 4.92                                 # wind speed [m/s]
    # F = 104956.99                               # dist where wind has const velocity [m/s] calculated from peak freq eq
    # alpha = 0.076 * (U_10**2 / (F * g))**0.22   # factor in spectral eq
    # w_p = (2*np.pi)/11                           # peak frequency [rad/s]
    # gamma = 1.54                                # enhancement factor

    # if w.any() <= w_p:
    #     sigma = 0.07
    # else:
    #     sigma = 0.09
    # r = np.exp(-(w - w_p)**2 / (2 * sigma**2 * w_p**2))
    # S_j = ((alpha*g**2)/w**5) * np.exp((-5/4) * (w_p/w)**4) * gamma**r 
    # print('S_j',S_j)