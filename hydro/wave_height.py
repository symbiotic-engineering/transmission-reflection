'''This script will find the wave elevation directly in front of and behind
each body. It will take the average of the total wave elevation in front of 
and behind each body and divide it by the average incident wave elevation
in front of and behind each body. This is how the reflection and transmission
coefficients are found, respectively. The energy blance and energy dissipation
are also computed.'''

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
    g = 9.81            # gravitational constant (m/s^2)
    k = w**2/g      # wave number infinite depth (rad^2/m)
    lam = int(2*np.pi/k)    # wavelength infinite depth (m)

    zinc_up = incoming_fse[mid_y, (mid_x - int(rel_dim*convx) + int(x_center*convx)) - int(lam*convx):mid_x - int(rel_dim*convx) + int(x_center*convx)]      # incident wave height upstream
    zinc_down = incoming_fse[mid_y, mid_x + int(rel_dim*convx) + int(x_center*convx):(mid_x + int(rel_dim*convx) + int(x_center*convx)) + int(lam*convx)]    # incident wave height downstream
    z_up = total[mid_y, (mid_x - int(rel_dim*convx) + int(x_center*convx)) - int(lam*convx):mid_x - int(rel_dim*convx) + int(x_center*convx)]                # total wave height upstream
    z_down = total[mid_y, mid_x + int(rel_dim*convx) + int(x_center*convx):(mid_x + int(rel_dim*convx) + int(x_center*convx)) + int(lam*convx)]              # transmitted wave height

    avg_H_zup = np.array([np.mean(abs(z_up))])
    avg_H_zincup = np.array([np.mean(abs(zinc_up))])
    avg_H_zdown = np.array([np.mean(abs(z_down))])
    avg_H_zincdown = np.array([np.mean(abs(zinc_down))])
    
    if farm:
        z_upWEC1 = total[mid_y + int(ytrans[0]*convy), (mid_x - int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx)) - int(lam*convx) :mid_x - int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx)]
        z_upWEC3 = total[mid_y + int(ytrans[1]*convy), (mid_x - int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx)) - int(lam*convx) :mid_x - int(rel_dim*convx) + int(xtrans[1]*convx) + int(x_center*convx)]

        z_downWEC1 = total[mid_y + int(ytrans[0]*convy), mid_x + int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx):(mid_x + int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx)) + int(lam*convx)]
        z_downWEC3 = total[mid_y + int(ytrans[1]*convy), mid_x + int(rel_dim*convx) + int(xtrans[1]*convx) + int(x_center*convx):(mid_x + int(rel_dim*convx) + int(xtrans[0]*convx) + int(x_center*convx)) + int(lam*convx)]

        avg_H_zup = np.array([np.mean(abs(z_upWEC1)), np.mean(abs(z_up)), np.mean(abs(z_upWEC3))])
        avg_H_zdown = np.array([np.mean(abs(z_downWEC1)), np.mean(abs(z_down)), np.mean(abs(z_downWEC3))])

    ref = abs(avg_H_zup / avg_H_zincup) - 1
    trans = abs(avg_H_zdown / avg_H_zincdown)
    EB = trans**2 + ref**2     # energy balance
    KD = 1 - EB  

    rho = 1025                  # [kg/m^3] density of sea water
    g = 9.81  
    # absorbed power per unit width
    power_abs = ((rho*g**2)/(4*w))*((avg_H_zincup/2)**2 - ((avg_H_zup - avg_H_zincup)/2)**2 - (avg_H_zdown/2)**2) # [W/m]
    print('absorber power',power_abs)

    return ref, trans, EB, KD, power_abs