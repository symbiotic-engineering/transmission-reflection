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

    if w.any() <= w_p:
        sigma = 0.07
    else:
        sigma = 0.09
    r = np.exp(-(w - w_p)**2 / (2 * sigma**2 * w_p**2))
    S_j = ((alpha*g**2)/w**5) * np.exp((-5/4) * (w_p/w)**4) * gamma**r 
    print('S_j',S_j)
    
    # S_j = []                                    # JONSWAP wave spectral density [m^2/Hz]

    # for omega in w:
    #     if omega <= w_p:
    #         sigma = 0.07
    #     else:
    #         sigma = 0.09
    #     r = np.exp(-(omega - w_p)**2 / (2 * sigma**2 * w_p**2))

    #     S = ((alpha*g**2)/omega**5) * np.exp((-5/4) * (w_p/omega)**4) * gamma**r 
    #     S_j.append(S)

    # f = w / (np.pi * 2)
    # plt.plot(f,S_j)
    # plt.xlabel('f [Hz]')
    # plt.ylabel('Spectral Density [m^2/Hz]')
    # plt.savefig('jonswap.pdf')

    return S_j

# # Main script
# w = np.linspace(0.5, 3.5, 100)  # Frequencies in rad/s
# f = w / (2 * np.pi)  # Convert to Hz
# gamma_values = np.array([3.3])  # Different gamma values

# # Adjust font sizes
# plt.rcParams.update({
#     'font.size': 18,  # Default font size
#     'axes.titlesize': 22,  # Title font size
#     'axes.labelsize': 20,  # Axes label font size
#     'xtick.labelsize': 18,  # X-tick label font size
#     'ytick.labelsize': 18,  # Y-tick label font size
#     'legend.fontsize': 18,  # Legend font size
#     'figure.titlesize': 22,  # Figure title font size
# })

# # Use a colorblind-friendly color palette
# colors = ["#0072B2", "#D55E00", "#CC79A7", "#009E73"]  # Blue, Vermilion, Yellow, Green

# plt.figure(figsize=(12, 8))
# for i, gamma in enumerate(gamma_values):
#     spectra = get_spectra(w, gamma)  # Compute spectra for all frequencies
#     plt.plot(w, spectra, label=f"γ = {gamma:.2f}", color=colors[i], linewidth=4)

# # Plot settings
# plt.xlabel('Frequency [rad/s]')
# plt.ylabel('Spectral Density [m^2/Hz]')
# plt.title('JONSWAP Spectra γ = 3.3')
# #plt.legend()
# plt.grid()
# plt.savefig('gtech_fig.pdf')
# plt.show()