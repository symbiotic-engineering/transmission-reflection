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

def farfield(point_absorber,oscillating_surge,controls,H,T,xgrid,ygrid,csv_file):
    mxc = int(xgrid/10)                                   # number of grid points in x (-1) (10 m res in x)
    myc = int(ygrid/10)                                   # number of grid points in y (-1) (10 m res in y)

    x1, x2, x3, x4 = 1490,1550,1460,1520
    #x = [x1,x2,x3,x4]
    x = [x1,x2,x3,x4,
         x1 + 200, x2 + 200, x3 + 200, x4 + 200,
         x1 - 200, x2 - 200, x3 - 200, x4 - 200,
         x1 + 400, x2 + 400, x3 + 400, x4 + 400,
         x1 - 400, x2 - 400, x3 - 400, x4 - 400]
    ya, yb = 4500, 4460

    ## to obtain Kt and Kr coefficients for your body and case
    ## note: won't run on my personal laptop, but SWAN won't run on lab computer
    # Kt_H, Kr_H, w_vals, power = run_coeffs.wec_run(w,breakwtr,point_absorber,oscillating_surge,attenuator,farm,controls,staggered,reactive)

    def configure_coefficients(csv_file):
        # Read the CSV file into a DataFrame
        df = pd.read_csv(csv_file)

        KT = [df.iloc[0, 0], df.iloc[0, 1], df.iloc[0, 2],df.iloc[0, 3],
              df.iloc[0, 0], df.iloc[0, 1], df.iloc[0, 2],df.iloc[0, 3],
              df.iloc[0, 0], df.iloc[0, 1], df.iloc[0, 2],df.iloc[0, 3],
              df.iloc[0, 0], df.iloc[0, 1], df.iloc[0, 2],df.iloc[0, 3],
              df.iloc[0, 0], df.iloc[0, 1], df.iloc[0, 2],df.iloc[0, 3]]
        KR = [df.iloc[0, 4], df.iloc[0, 5], df.iloc[0, 6],df.iloc[0, 7],
              df.iloc[0, 4], df.iloc[0, 5], df.iloc[0, 6],df.iloc[0, 7],
              df.iloc[0, 4], df.iloc[0, 5], df.iloc[0, 6],df.iloc[0, 7],
              df.iloc[0, 4], df.iloc[0, 5], df.iloc[0, 6],df.iloc[0, 7],
              df.iloc[0, 4], df.iloc[0, 5], df.iloc[0, 6],df.iloc[0, 7]]

        # Print results
        print("KR:", KR)
        print("KT:", KT)

        # Return the KR and KT arrays
        return KR, KT

    KR, KT = configure_coefficients(csv_file)

    # diameter of the bodies
    if point_absorber:
        d = 14.6
    if oscillating_surge:
        d = 18

    # run SWAN and generate wave height data
    sfgrid_dat, sfgrid_tbl = run_swan.generate_swan_input(KR, KT, d, x, ya, yb, H, T, xgrid, ygrid, mxc, myc)

    return sfgrid_dat