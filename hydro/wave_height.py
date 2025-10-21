'''This script will find the wave elevation directly in front of and behind
each body. It will take the average of the total wave elevation in front of 
and behind each body and divide it by the average incident wave elevation
in front of and behind each body. This is how the reflection and transmission
coefficients are found, respectively. The energy blance, energy dissipation,
and power per unit width are also computed.'''

def wave_height(total,incoming_fse,xtrans,ytrans,rel_dim,w,nx,ny,x1,x2,y1,y2,x_center):
    import numpy as np
    import matplotlib.pyplot as plt
    import warnings
    warnings.filterwarnings("ignore", category=np.ComplexWarning)
    
    ##################################################################
    # Extract z_up and z_down using rel_dim to determine the region
    #                               --> need to account for x-shift in body position
    #                               --> need to ensure still averaging over a wavelength
    g = 9.81                # gravitational constant (m/s^2)
    k = w**2/g              # wave number infinite depth (rad^2/m)
    rho = 1025                                  # [kg/m^3] density of sea water
        
    # x- and y-positions, converted to grid points
    convx = int(nx/(np.abs(x1)+x2))      # to convert meters to grid points
    convy = int(ny/(abs(y1)+y2))

    x_start, y_start = x_center, 0
    Dx, Dy = xtrans, ytrans

    posx = np.array([x_start,x_start,x_start + Dx,x_start + Dx])*convx 
    posy = np.array([y_start, y_start + Dy,y_start + (1/2)*Dy,y_start + (3/2)*Dy])*convy

    # other misc conversions
    rel_dim = rel_dim*convx
    lam = int(2*np.pi/k)*convx    # wavelength infinite depth (m)

    zinc_up, zinc_down, z_up, z_down = [],[],[],[]
    ref, trans, EB, KD = [],[],[],[]
    power_abs = []

    for i in range(4):
        zinc_up_loop = incoming_fse[int(posy[i]), int(posx[i] - rel_dim - lam) : int(posx[i] - rel_dim)]      # incident wave height upstream
        zinc_down_loop = incoming_fse[int(posy[i]), int(posx[i] + rel_dim) : int(posx[i] + rel_dim + lam)]    # incident wave height downstream
        z_up_loop = total[int(posy[i]), int(posx[i] - rel_dim - lam) : int(posx[i] - rel_dim)]                # total wave height upstream
        z_down_loop = total[int(posy[i]), int(posx[i] + rel_dim) : int(posx[i] + rel_dim + lam)]              # transmitted wave height

        zinc_up.append(zinc_up_loop)
        zinc_down.append(zinc_down_loop)
        z_up.append(z_up_loop)
        z_down.append(z_down_loop)

        avg_incUP = np.mean(np.abs(zinc_up_loop))
        avg_incDOWN = np.mean(np.abs(zinc_down_loop))
        avg_ref = np.mean((np.abs(z_up_loop) - np.abs(zinc_up_loop)))
        avg_trans = np.mean(np.abs(z_down_loop))

        ref_loop = (avg_ref / avg_incUP)                 # reflection coefficient
        trans_loop = avg_trans / avg_incDOWN             # transmission coefficient
        EB_loop = trans_loop**2 + ref_loop**2                      # energy balance
        KD_loop = 1 - EB_loop                                 # dissipation coefficient

        ref.append(ref_loop)
        trans.append(trans_loop)
        EB.append(EB_loop)
        KD.append(KD_loop)
    
        # absorbed power per unit width
        power_absloop = ((rho*g**2)/(4*w))*((avg_incUP/2)**2 - ((avg_ref)/2)**2 - (avg_trans/2)**2) # [W/m]
        power_abs.append(power_absloop)

    print('absorbed power',power_abs)

    upstreamx = np.linspace(int(posx[1] - rel_dim - lam)/convx,int(posx[1] - rel_dim)/convx,num = np.size(zinc_up[1]))
    downstreamx = np.linspace(int(posx[1] + rel_dim)/convx,int(posx[1] + rel_dim + lam)/convx,num=np.size(zinc_down[1]))

    plt.plot(upstreamx, np.abs(zinc_up[0]),label='zinc_up',marker='o')
    plt.plot(downstreamx, np.abs(zinc_down[0]),label='zinc_down',marker='x')
    plt.plot(upstreamx, np.abs(z_up[0]),label='z_up',marker='s')
    plt.plot(downstreamx, np.abs(z_down[0]),label='z_down',marker='p')
    plt.legend()
    print('tic')
    plt.savefig('wave_heights.pdf')
    print('toc')
    plt.clf()

    return ref, trans, EB, KD, power_abs