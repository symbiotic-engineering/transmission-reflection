'''This is the run file to run all the near field functions. It runs:
1. body.name() which initializes the body
2. solve.hydro() which solves the body hydrodynamics (and includes the
power take-off function nested within)
3. solve.elevation() which solves for the free surface elevation
4. wave_height.wave_height() which solves for wave height ratios i.e.
the transmission and reflection coefficients.'''

'''LAST UPDATE" OCT 16TH 2025'''

def wec_run(w,point_absorber,oscillating_surge,controls):
    ################ SCRIPT FOR CALCULATING AT COEFFS ############
    # run all Kt and Kr calcs for any body here
    import body         # "body" contains functions for PA and OSWEC
    import solve        # this solves the hydrodynamics
    import wave_height  # this function finds Kt and Kr based on wave elevation
    import numpy as np
    import matplotlib.pyplot as plt 

    index = 4
    B = 0                                           # wave direction [rad]
    depth = 500                                     # keep deep water assumption for EB
    scale = 50

    xtrans = 40
    ytrans = 40                     # y translation of bodies
    x_center = 300                  # arbitrary x_center position to make indexing easier later

    Kr_H = [[] for _ in range(index)]       # initializing reflection coeff
    Kt_H = [[] for _ in range(index)]       # initializing transmission coeff
    power = [[] for _ in range(index)]      # initializing absorberd power [W/m]

    w_vals = []                             # for storing omega values

    if point_absorber:
        BdPA_num = 13.2332*w**2 + 10.7815*w + 1.8630 
        BdPA_den = w**2 + 1.0961*w + 0.3752
        # empirical damping fit for isolated device
        B_diff = (scale**(5/2))*(BdPA_num/BdPA_den) * np.exp(1j*(-0.0176*w**2 + 0.1140*w + 1.3173))

        # empirical damping fit for damping introduced to each device in the array
        Bd1 = (scale**(5/2))*(270.8966 + (-15.7783 - 270.8966)/(1 + (w/5.3854)**38.6287)) + 1j * (634.7380 + (263.8006 - 634.7380)/(1 + (w/5.0931)**108.8870))
        Bd2 = (scale**(5/2))*(25.0043 + (-16.6685 - 25.0043)/(1 + (w/5.4484)**71.1566)) + 1j * (91.9514 + 25.3225*np.cos(w*2.7611) + 9.7855*np.sin(w*2.7611))
        Bd3 = (scale**(5/2))*(-70.2006 - 51.1202*np.cos(w*1.4738) - 22.1776*np.sin(w*1.4738)) + 1j * (-126.9010 - 6.5431*np.cos(w*3.0603) + 19.0480*np.cos(w*3.0603))
        Bd4 = (scale**(5/2))*(811.5504*w**3 - 1.2703e+04*w**2 + 6.5620e+04*w - 1.1218e+05) + 1j * (-110.5101 + (-507.6379 + 110.5101)/(1 + (w/5.2618)**-57.2704))
    else:
        # empirical damping fit for isolated device
        B_diff = (scale**(5/2))*(0.7731*w**2 - 6.1100*w + 13.3851) * np.exp(1j*(0.0219*w**2 - 0.19911*w + 0.6697))

        # empirical damping fit for damping introduced to each device in the array
        Bd1 = (scale**(5/2))*(22.3463*w**3 - 361.6316*w**2 + 1.9487e+03*w - 3.4691e+03) + 1j * (12.2825*w**2 - 95.4165*w + 245.6894)
        Bd2 = (scale**(5/2))*(8.4637 + 57.0776*np.cos(w*3.9128) + 21.1357*np.sin(w*3.9128)) + 1j * (21.7659*w**2 - 193.2818*w + 491.0362)
        Bd3 = (scale**(5/2))*(3.7682*w**2 - 69.7304*w + 229.6366) + 1j * (17.8272*w**3 - 309.0917*w**2 + 1.7637e+03*w - 3.2040e+03)
        Bd4 = (scale**(5/2))*(-45.0017*w**3 + 719.5501*w**2 - 3.8389e+03*w + 6.7988e+03) + 1j * (95.3911*w**3 - 1.5521e+03*w**2 + 8.3492e+03*w - 1.4736e+04)
    
    for i in range(np.size(w)):
        res = 1                                                                                                  # set the grid resolution
        array, rel_dim, char_dim = body.initialize(xtrans,ytrans,w[i],x_center,point_absorber)                                      # generate the meshed array
        diff_result,rad_result,RAO_vals,lam,CWR = solve.hydro(array,B,depth,w[i],char_dim,controls,point_absorber,B_diff,Bd1[i],Bd2[i],Bd3[i],Bd4[i])  # solve hydrodynamics
        total,incoming_fse,x1,x2,nx,y1,y2,ny = solve.elevation(res,lam,diff_result,rad_result,RAO_vals,controls,rel_dim) # solve for wave elevation
        ref,trans,EB,KD,power_abs = wave_height.wave_height(total,incoming_fse,xtrans,ytrans,rel_dim,w,nx,ny,x1,x2,y1,y2,x_center) # calculate reflection and transmission coefficients

        print('Kt',trans)
        print('Kr',ref)
        print('dissipation',KD)

        w_vals.append(w)

        for i in range(index):
            Kr_H[i].append(ref[i])
            Kt_H[i].append(trans[i])
            power[i].append(power_abs[i])
        
    return Kt_H, Kr_H, w_vals, power, RAO_vals