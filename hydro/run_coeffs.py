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
    import body         # "body" contains functions for PA, OSWEC, attenuator, and breakwater initialization
    import solve        # this solves the hydrodynamics
    import wave_height  # this function finds Kt and Kr based on wave elevation
    import numpy as np
    import matplotlib.pyplot as plt 

    N = 3
    index = 3

    B = 0                                           # wave direction [rad]
    depth = 500                                     # keep deep water assumption for EB
    
    if staggered:
        xtrans = np.array([50,50])                  # x translation of bodies if staggered farm
        x_center = -25
    else:
        xtrans = np.array([0,0])
        x_center = 0

    ytrans = np.array([50,-50])                     # y translation of bodies if farm

    Kr_H = [[] for _ in range(index)]       # initializing reflection coeff
    Kt_H = [[] for _ in range(index)]       # initializing transmission coeff
    power = [[] for _ in range(index)]      # initializing absorberd power [W/m]

    w_vals = []                             # for storing omega values

    for w in w:
        # this is where you set the grid resolution based on wavelength.
        # shorter wavelengths require finer mesh resolution based on 
        # mesh convergence study.
        if w < 1.0:
            res = 2.0
        else:
            res = 3.5

        # this where the code generates the body based on which you chose
        if point_absorber:
            array, rel_dim, char_dim = body.PA(xtrans,ytrans,w,x_center)
        if oscillating_surge:
            array, rel_dim, char_dim = body.OSWEC(xtrans,ytrans,w,x_center)

        # this is where the code solves hydrodynamics
        diff_result,rad_result,RAO_vals,lam,CWR = solve.hydro(array,B,depth,w,char_dim,controls,point_absorber)

        # this is where the code solves for wave elevation
        total,incoming_fse,x1,x2,nx,y1,y2,ny = solve.elevation(res,lam,diff_result,rad_result,RAO_vals,controls,N,rel_dim)

        # this is where the code calculates your reflection and transmission coefficients
        ref,trans,EB,KD,power_abs = wave_height.wave_height(total,incoming_fse,xtrans,ytrans,rel_dim,w,nx,ny,x1,x2,y1,y2,x_center)

        print('Kt',trans)
        print('Kr',ref)
        print('dissipation',KD)

        w_vals.append(w)

        for i in range(index):
            Kr_H[i].append(ref[i])
            Kt_H[i].append(trans[i])
            power[i].append(power_abs[i])
        
    return Kt_H, Kr_H, w_vals, power, RAO_vals