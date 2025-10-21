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

    N = 3
    index = 4

    B = 0                                           # wave direction [rad]
    depth = 500                                     # keep deep water assumption for EB

    xtrans = 40
    ytrans = 40                     # y translation of bodies
    x_center = 300                  # arbitrary x_center position to make indexing easier later

    Kr_H = [[] for _ in range(index)]       # initializing reflection coeff
    Kt_H = [[] for _ in range(index)]       # initializing transmission coeff
    power = [[] for _ in range(index)]      # initializing absorberd power [W/m]

    w_vals = []                             # for storing omega values

    for w in w:

        res = 1                                                                                                  # set the grid resolution
        array, rel_dim, char_dim = body.initialize(xtrans,ytrans,w,x_center,point_absorber)                                      # generate the meshed array
        diff_result,rad_result,RAO_vals,lam,CWR = solve.hydro(array,B,depth,w,char_dim,controls,point_absorber)  # solve hydrodynamics
        total,incoming_fse,x1,x2,nx,y1,y2,ny = solve.elevation(res,lam,diff_result,rad_result,RAO_vals,controls,N,rel_dim) # solve for wave elevation
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