'''This script will find the wave elevation directly in front of and behind
each body. It will take the average of the total wave elevation in front of 
and behind each body and divide it by the average incident wave elevation
in front of and behind each body. This is how the reflection and transmission
coefficients are found, respectively. The energy blance, energy dissipation,
and power per unit width are also computed.'''

def wave_height(total, incoming_fse, xtrans, ytrans, farm, rel_dim, w, nx, ny, x1, x2, y1, y2,x_center):
    import numpy as np
    import matplotlib.pyplot as plt
    import warnings
    warnings.filterwarnings("ignore", category=np.ComplexWarning)
    
    ##################################################################
    # Extract relevant columns
    mid_y = int(ny / 2)
    mid_x = int(nx / 2)
    convx = int(nx/(abs(x1)+x2))      # to convert meters to grid points
    convy = int(ny/(abs(y1)+y2))
    
    # Extract z_up and z_down using rel_dim to determine the region
    #                               --> need to account for x-shift in body position
    #                               --> need to ensure still averaging over a wavelength
    g = 9.81                # gravitational constant (m/s^2)
    k = w**2/g              # wave number infinite depth (rad^2/m)
    lam = int(2*np.pi/k)    # wavelength infinite depth (m)

    zinc_up = incoming_fse[mid_y, (mid_x - int(rel_dim*convx) + int(x_center*convx)) - int(lam*convx):mid_x - int(rel_dim*convx) + int(x_center*convx)]      # incident wave height upstream
    zinc_down = incoming_fse[mid_y, mid_x + int(rel_dim*convx) + int(x_center*convx):(mid_x + int(rel_dim*convx) + int(x_center*convx)) + int(lam*convx)]    # incident wave height downstream
    z_up = total[mid_y, (mid_x - int(rel_dim*convx) + int(x_center*convx)) - int(lam*convx):mid_x - int(rel_dim*convx) + int(x_center*convx)]                # total wave height upstream
    z_down = total[mid_y, mid_x + int(rel_dim*convx) + int(x_center*convx):(mid_x + int(rel_dim*convx) + int(x_center*convx)) + int(lam*convx)]              # transmitted wave height

    avg_H_zup = np.array([np.mean(abs(z_up - zinc_up))])
    avg_H_zincup = np.array([np.mean(abs(zinc_up))])
    avg_H_zdown = np.array([np.mean(abs(z_down))])
    avg_H_zincdown = np.array([np.mean(abs(zinc_down))])
    
    if farm:
        zinc_upWEC1 = incoming_fse[mid_y + int(ytrans[0]*convy), (mid_x - int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx)) - int(lam*convx) :mid_x - int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx)]
        zinc_upWEC3 = incoming_fse[mid_y + int(ytrans[1]*convy), (mid_x - int(rel_dim*convx) + int(xtrans[1]*convx) + int(x_center*convx)) - int(lam*convx) :mid_x - int(rel_dim*convx) + int(xtrans[1]*convx) + int(x_center*convx)]
        z_upWEC1 = total[mid_y + int(ytrans[0]*convy), (mid_x - int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx)) - int(lam*convx) :mid_x - int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx)]
        z_upWEC3 = total[mid_y + int(ytrans[1]*convy), (mid_x - int(rel_dim*convx) + int(xtrans[1]*convx) + int(x_center*convx)) - int(lam*convx) :mid_x - int(rel_dim*convx) + int(xtrans[1]*convx) + int(x_center*convx)]

        zinc_downWEC1 = incoming_fse[mid_y + int(ytrans[0]*convy), mid_x + int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx):(mid_x + int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx)) + int(lam*convx)]
        zinc_downWEC3 = incoming_fse[mid_y + int(ytrans[1]*convy), mid_x + int(rel_dim*convx) + int(xtrans[1]*convx) + int(x_center*convx):(mid_x + int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx)) + int(lam*convx)]
        z_downWEC1 = total[mid_y + int(ytrans[0]*convy), mid_x + int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx):(mid_x + int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx)) + int(lam*convx)]
        z_downWEC3 = total[mid_y + int(ytrans[1]*convy), mid_x + int(rel_dim*convx) + int(xtrans[1]*convx) + int(x_center*convx):(mid_x + int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx)) + int(lam*convx)]

        avg_H_zincup = np.array([np.mean(abs(zinc_upWEC1)), np.mean(abs(zinc_up)), np.mean(abs(zinc_upWEC3))])
        avg_H_zincdown = np.array([np.mean(abs(zinc_downWEC1)), np.mean(abs(zinc_down)), np.mean(abs(zinc_downWEC3))])
        avg_H_zup = np.array([np.mean(abs(z_upWEC1 - zinc_upWEC1)), np.mean(abs(z_up - zinc_up)), np.mean(abs(z_upWEC3 - zinc_upWEC3))])
        avg_H_zdown = np.array([np.mean(abs(z_downWEC1)), np.mean(abs(z_down)), np.mean(abs(z_downWEC3))])

    ### reflection coeff note: you get the same value if you subtract complex incident
    ### wave elevation from complex upstream wave, then dividing abs values of new
    ### upstream wave by abs value of incident wave as when you
    ### divide abs value of total upstream wave by abs value of incident wave
    ### and then subtract one

    ref = (avg_H_zup / avg_H_zincup) # - 1
    trans = avg_H_zdown / avg_H_zincdown
    EB = trans**2 + ref**2     # energy balance
    KD = 1 - EB  

    rho = 1025                  # [kg/m^3] density of sea water
    g = 9.81  
    # absorbed power per unit width
    power_abs = ((rho*g**2)/(4*w))*((avg_H_zincup/2)**2 - ((avg_H_zup)/2)**2 - (avg_H_zdown/2)**2) # [W/m]
    print('absorber power',power_abs)

    return ref, trans, EB, KD, power_abs