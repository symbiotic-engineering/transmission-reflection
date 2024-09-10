'''this script will run "run_coeffs" for all your specified cases 
and save the .csv files in the "data" folder. This data can then 
be easily used for post processing in the "post_pro" folder'''

import numpy as np
import matplotlib.pyplot as plt
import csv
import run_coeffs
import os

w = np.array([0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.3])  # wave frequency

# Define the combinations of True and False for your parameters
combinations = [
    {'breakwtr': False, 'point_absorber': True, 'oscillating_surge': False, 'attenuator': False, 'farm': True, 'controls': False, 'staggered': False, 'reactive': False},
    {'breakwtr': False, 'point_absorber': True, 'oscillating_surge': False, 'attenuator': False, 'farm': True, 'controls': True, 'staggered': False, 'reactive': False},
    {'breakwtr': False, 'point_absorber': True, 'oscillating_surge': False, 'attenuator': False, 'farm': True, 'controls': True, 'staggered': False, 'reactive': True},
    {'breakwtr': False, 'point_absorber': True, 'oscillating_surge': False, 'attenuator': False, 'farm': True, 'controls': False, 'staggered': True, 'reactive': False},
    {'breakwtr': False, 'point_absorber': True, 'oscillating_surge': False, 'attenuator': False, 'farm': True, 'controls': True, 'staggered': True, 'reactive': False},
    {'breakwtr': False, 'point_absorber': True, 'oscillating_surge': False, 'attenuator': False, 'farm': True, 'controls': True, 'staggered': True, 'reactive': True}
]

# Map combinations to their corresponding file names
# breakwtr,point_absorber,oscillating_surge,attenuator,farm,controls,staggered,reactive
file_names = {
    (False, True, False, False, True, False, False, False): 'PA_reg_uncont.csv',
    (False, True, False, False, True, True, False, False): 'PA_reg_damp.csv',
    (False, True, False, False, True, True, False, True): 'PA_reg_react.csv',
    (False, True, False, False, True, False, True, False): 'PA_stag_uncont.csv',
    (False, True, False, False, True, True, True, False): 'PA_stag_damp.csv',
    (False, True, False, False, True, True, True, True): 'PA_stag_react.csv'
}

# Loop through each combination
for combination in combinations:
    # Unpack the dictionary to pass the parameters to your function
    Kt_H, Kr_H, w_vals, power = run_coeffs.wec_run(w, **combination)

    # Create a tuple of the current combination to use as a key for the file name
    combination_key = tuple(combination.values())

    data_folder = os.path.join(os.path.dirname(__file__), 'data')
    file_path = os.path.join(data_folder, file_names.get(combination_key, 'default.csv'))

    # Save the data to the corresponding .csv file
    with open(file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        header = ['Omega'] + [f'Kt_H_{i+1}' for i in range(len(Kt_H))] + [f'Kr_H_{i+1}' for i in range(len(Kr_H))] + [f'power_{i+1}' for i in range(len(power))]
        writer.writerow(header)
        for i in range(len(w_vals)):
            row = [w_vals[i]] + [Kt_H[j][i] for j in range(len(Kt_H))] + [Kr_H[j][i] for j in range(len(Kr_H))] + [power[j][i] for j in range(len(power))]
            writer.writerow(row)
    
    print('RUN COMPLETE:',file_path)
