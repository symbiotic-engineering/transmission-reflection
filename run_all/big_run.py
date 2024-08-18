''' TO RUN THIS SCRIPT (read the README too): 
    1. You need to define your xgrid and ygrid size and resolution.
    2. You need to check the xtrans and ytrans WEC positions in the
    run_coeffs.wec_run() function. This will determine if you will run a 
    staggered or regular array.
    3. You need to define the positions of your bodies in the SWAN grid
    according to the run_coeffs.wec_run() function.
    4. You need to define the significant wave height and peak wave
    period of your sea state. ATTN: you will also need to adjust the
    run_swan.py file for wind speed and direction.
    5. Decide which body you are running by setting it equal to True.
    Set all others to False.
    6. If you want to include array interactions, set farm equal to
    True. If you want to include controls, set controls equal to True.
    If you want reactive controls, set reactive equal to True; otherwise,
    if controls == True and reactive == False, it will do damped control.'''

import sys
import os
import numpy as np
import pandas as pd
# Get the current directory of the script
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
hydro_dir = os.path.join(parent_dir, 'hydro')
swan_dir = os.path.join(parent_dir, 'swan')
sys.path.append(hydro_dir)
sys.path.append(swan_dir)
import run_coeffs
import gen_data
import run_swan

def farfield(six,swan_stag,breakwtr,point_absorber,oscillating_surge,attenuator,farm,controls,staggered,reactive,H,T,xgrid,ygrid,csv_file):
    mxc = int(xgrid/10)                                   # number of grid points in x (-1) (10 m res in x)
    myc = int(ygrid/10)                                   # number of grid points in y (-1) (10 m res in y)

    if six:
        if swan_stag:
            # six body stag: 
            x = [1490,1550,1610, 1440, 1500, 1560]
            ya, yb = 4500, 4450
        else:
            # six-body reg: 
            x = [1290,1350,1410,1470,1530,1590]    # x-position of bodies
            ya, yb = 4500,4500                     # y-position of rows 1 and 2
    else:
        if swan_stag:
            # three-body stag: 
            x = [1490,1550,1610] 
            ya, yb = 4500, 4450         # change i==1 in run_swan
        else:
            # three-body reg: 
            x = [1490,1550,1610] 
            ya, yb = 4500, 4450

    ## to obtain Kt and Kr coefficients for your body and case
    ## note: won't run on my personal laptop, but SWAN won't run on lab computer
    # Kt_H, Kr_H, w_vals, power = run_coeffs.wec_run(w,breakwtr,point_absorber,oscillating_surge,attenuator,farm,controls,staggered,reactive)

    def configure_coefficients(csv_file, farm):
        # Read the CSV file into a DataFrame
        df = pd.read_csv(csv_file)
        if farm:
            i = 1
        else:
            i = 0

        if farm:
            Kt_H = [df.iloc[5, 1], df.iloc[5, 2], df.iloc[5, 3]]
            Kr_H = [df.iloc[5, 4], df.iloc[5, 5], df.iloc[5, 6]]
        else:
            Kt_H = [df.iloc[5, 1]]
            Kr_H = [df.iloc[5, 2]]

        # Configure the KR and KT arrays
        KR = [Kr_H[0], Kr_H[i], Kr_H[2*i], Kr_H[0], Kr_H[i], Kr_H[2*i]]
        KT = [Kt_H[0], Kt_H[i], Kt_H[2*i], Kt_H[0], Kt_H[i], Kt_H[2*i]]

        # Print results
        print("KR:", KR)
        print("KT:", KT)

        # Return the KR and KT arrays
        return KR, KT

    KR, KT = configure_coefficients(csv_file, farm)

    # diameter of the bodies
    if point_absorber:
        d = 21
    if oscillating_surge:
        d = 18
    if breakwtr:
        d = 20
    if attenuator:
        d = 4

    # run SWAN and generate wave height data
    sfgrid_dat, sfgrid_tbl = run_swan.generate_swan_input(KR, KT, d, x, ya, yb, H, T, xgrid, ygrid, mxc, myc,attenuator,six,swan_stag)

    return sfgrid_dat