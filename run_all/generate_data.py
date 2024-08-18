# Define the combinations and their corresponding file names
# 1: body's coeffs found in isolation                           a: three-body reg
# 2: body's coeffs found in reg array                           b: three-body stag
# 3: body's coeffs found in stag array                          c: six-body reg
#                                                               d: six-body stag

# six, swan_stag, breakwtr, point_absorber, oscillating_surge, attenuator, farm, controls, staggered, reactive,

import numpy as np
import sys
import os
import pandas as pd
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
swan_dir = os.path.join(parent_dir, 'swan')
sys.path.append(swan_dir)
import gen_data
import big_run

# Define the combinations and their corresponding file names
file_map = {
    # breakwater
    # (False, False, True, False, False, False, False, False, False, False): ('break_1a.csv', 'break.csv'),
    # (False, True, True, False, False, False, False, False, False, False): ('break_1b.csv', 'break.csv'),
    # (True, False, True, False, False, False, False, False, False, False): ('break_1c.csv', 'break.csv'),
    # (True, True, True, False, False, False, False, False, False, False): ('break_1d.csv', 'break.csv'),
    # (False, False, True, False, False, False, True, False, False, False): ('break_2a.csv', 'break_reg.csv'),
    # (False, True, True, False, False, False, True, False, True, False): ('break_3b.csv', 'break_stag.csv')
    # point absorber
    # (False, False, False, True, False, False, False, True, False, False): ('PA_1a.csv', 'PA_damp.csv'),
    # (False, True, False, True, False, False, False, True, False, False): ('PA_1b.csv', 'PA_damp.csv'),
    # (True, False, False, True, False, False, False, True, False, False): ('PA_1c.csv', 'PA_damp.csv'),
    # (True, True, False, True, False, False, False, True, False, False): ('PA_1d.csv', 'PA_damp.csv'),
    # (False, False, False, True, False, False, True, True, False, False): ('PA_2a.csv', 'PA_reg_damp.csv'),
    # (False, True, False, True, False, False, True, True, True, False): ('PA_3b.csv', 'PA_stag_damp.csv'),
    # # oscillating surge
    # (False, False, False, False, True, False, False, True, False, False): ('OS_1a.csv', 'OS_damp.csv'),
    # (False, True, False, False, True, False, False, True, False, False): ('OS_1b.csv', 'OS_damp.csv'),
    # (True, False, False, False, True, False, False, True, False, False): ('OS_1c.csv', 'OS_damp.csv'),
    # (True, True, False, False, True, False, False, True, False, False): ('OS_1d.csv', 'OS_damp.csv'),
    # (False, False, False, False, True, False, True, True, False, False): ('OS_2a.csv', 'OS_reg_damp.csv'),
    # (False, True, False, False, True, False, True, True, True, False): ('OS_3b.csv', 'OS_stag_damp.csv'),
    # # attenuator
    # (False, False, False, False, False, True, False, True, False, False): ('atten_1a.csv', 'atten_damp.csv'),
    # (False, True, False, False, False, True, False, True, False, False): ('atten_1b.csv', 'atten_damp.csv'),
    # (True, False, False, False, False, True, False, True, False, False): ('atten_1c.csv', 'atten_damp.csv'),
    # (True, True, False, False, False, True, False, True, False, False): ('atten_1d.csv', 'atten_damp.csv'),
    # (False, False, False, False, False, True, True, True, False, False): ('atten_2a.csv', 'atten_reg_damp.csv'),
    # (False, True, False, False, False, True, True, True, True, False): ('atten_3b.csv', 'atten_stag_damp.csv')
    # blank dataset
    (False, True, False, False, False, True, True, True, True, False): ('blank.csv', 'atten_stag_damp.csv')
}

# input parameters from SouthFork Wind location                           
H = 0.8                                     # avg significant wave height [m]
T = 5                                       # avg wave period [s] from buoy 44097

# define computational grid
xgrid = 3000                                # size of grid in x-direction
ygrid = 5000                                # size of grid in y-direction

# Define the data directory where the files will be saved
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # transmission_reflection
data_dir = os.path.join(base_dir, 'data')

# Loop through each combination in file_map and run the cases
for key, (output_filename, csv_file_name) in file_map.items():
    six, swan_stag, breakwtr, point_absorber, oscillating_surge, attenuator, farm, controls, staggered, reactive = key

    # Construct the path to the output .csv file in transmission_reflection/data/
    output_filepath = os.path.join(data_dir, output_filename)

    # Construct the path to the CSV file being read (in hydro/data/)
    hydro_data_dir = os.path.join(base_dir, 'hydro', 'data')
    csv_file = os.path.join(hydro_data_dir, csv_file_name)

    sfgrid_dat = big_run.farfield(six, swan_stag, breakwtr, point_absorber, oscillating_surge, attenuator, farm, controls, staggered, reactive, H, T, xgrid, ygrid, csv_file)

    # Save .csv file in transmission_reflection/data/ and plot wave height data
    h_s = gen_data.save_wave_height_to_csv(sfgrid_dat, output_filepath)

    # Optional: Add any additional processing or logging here


#waveHeight = post_process.postpro(sfgrid_dat,xgrid,ygrid,mxc,myc,x,ya,yb)