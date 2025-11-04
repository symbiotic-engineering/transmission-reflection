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
    # point absorber
    (True,False,True): ('PA_SF.csv', 'PA_SF_damp.csv'),
    #(True,False,True): ('PA_large.csv','PA_damp.csv')
    # oscillating surge
    (False, True, True): ('OS_SF.csv', 'OS_SF_damp.csv')
    #(False, False, False, False, True, False, True, True, False, False): ('OS_2.csv', 'OS_uncont.csv')
}

# input parameters from SouthFork Wind location                           
H = 2.19                                        # avg significant wave height [m]
T = 11.0                                        # dominant wave period [s] from buoy 44097

# define computational grid
xgrid = 3000                                # size of grid in x-direction
ygrid = 5000                                # size of grid in y-direction

# Define the data directory where the files will be saved
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # transmission_reflection
data_dir = os.path.join(base_dir, 'data')

# Loop through each combination in file_map and run the cases
for key, (output_filename, csv_file_name) in file_map.items():
    point_absorber,oscillating_surge,controls = key

    # Construct the path to the output .csv file in transmission_reflection/data/
    output_filepath = os.path.join(data_dir, output_filename)

    # Construct the path to the CSV file being read (in hydro/data/)
    hydro_data_dir = os.path.join(base_dir, 'hydro', 'data')
    csv_file = os.path.join(hydro_data_dir, csv_file_name)

    sfgrid_dat = big_run.farfield(point_absorber, oscillating_surge,controls,H, T, xgrid, ygrid, csv_file)

    # Save .csv file in transmission_reflection/data/ and plot wave height data
    h_s = gen_data.save_wave_height_to_csv(sfgrid_dat, output_filepath)

    # Optional: Add any additional processing or logging here


#waveHeight = post_process.postpro(sfgrid_dat,xgrid,ygrid,mxc,myc,x,ya,yb)