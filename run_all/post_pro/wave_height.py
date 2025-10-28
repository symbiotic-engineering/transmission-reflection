import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

# Constants
mxc = 300
myc = 500
xgrid = 3000  # Total width in meters
ygrid = 5000  # Total height in meters

# Load and reshape data
def load_and_reshape(file_path):
    data = np.loadtxt(file_path, delimiter=',')
    return data.reshape(myc + 1, mxc + 1)

# Set the correct data folder path
data_folder = '../../data/'

# Define file paths
file_paths = [
    'blank.csv',
    #'OSr_1a.csv', 
    #'OSr_2a.csv',
    'PA_large.csv'
    #'PA_2.csv'
]

# Load data arrays
data_arrays = {fp: load_and_reshape(os.path.join(data_folder, fp)) for fp in file_paths}

# Calculate percent difference compared to blank.csv
blank_data = data_arrays['blank.csv']
percent_differences = {}

for fp in file_paths[1:]:  # Skip the 'blank.csv'
    percent_differences[fp] = abs(100 * (data_arrays[fp] - blank_data) / blank_data)

# Prepare data for saving
results_list = []

# Define the nautical miles of interest and corresponding indices
nautical_miles = np.array([0.25, 0.5, 1.0, 1.5, 2])
nautical_miles_conversion = 1852  # Conversion factor from meters to nautical miles

# for plotting only
# Define colorblind-friendly colors and linestyles
cud_colors = ['#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7', '#000000','#F5C200']
linestyles = ['-', '--', '-.', ':', '-', '--', '-.', ':']  # Make sure we have enough linestyles

#[1490,1550,1460,1520]
# Plot percent difference at x_center location from y_start to y=0 for all files
plt.figure()        #figsize=(8, 8))
for i, file_path in enumerate(file_paths[1:]):  # Skip the 'blank.csv'
    x_center = int((1460+1550)/2 * (mxc / xgrid))
    y_start = int(4450 * (myc / ygrid))
    
    # Extract percent difference at x_center from y_start to y = 0
    x_range = np.linspace(0,xgrid,mxc)*mxc/xgrid
    nautical_string = ['2.0 nm','1.5 nm','1.0 nm','0.5 nm','0.25 nm']

        #y_range = np.linspace(y_start * (ygrid / myc), 0, y_start + 1)
    plt.figure(figsize=(8,6))
    for j in range(np.size(nautical_miles)):
        y_dists_nautical = int(nautical_miles[j] * nautical_miles_conversion * myc/ygrid)
        percent_diff = percent_differences[file_path][y_dists_nautical, int(x_range[0]):int(np.size(x_range))]
        
        # Reverse the order for the plot
        # y_range = y_range[::-1] / 1852  # Convert to nautical miles
        # percent_diff = percent_diff[::-1]

        #parts = file_path.split('_')
        label = nautical_string[j]            #f"{parts[1].replace('.csv', '')}"
        
        # Plot the percent difference with colorblind-friendly colors and different linestyles
        plt.plot(x_range*xgrid/mxc, percent_diff, label=label,
                color=cud_colors[j])             #, linestyle=linestyles[i % len(linestyles)])

# Customize the plot
plt.xlabel('x [m]',fontsize=20)
#plt.ylabel('Wave Height Reduction [%]',fontsize=20)
plt.legend(loc='upper right',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylim([0,10])
ax = plt.gca()
ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
#plt.text(1.025, 1.05, 'a', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top', ha='left')
plt.tight_layout()
plt.savefig('PAlarge_reduction.pdf')
plt.show()

#     # Extract percent difference at x_center from y_start to y = 0
#     y_range = np.linspace(y_start * (ygrid / myc), 0, y_start + 1)
#     percent_diff_at_x_center = percent_differences[file_path][0:y_start + 1, x_center]
    
#        # Reverse the order for correct alignment
#     y_range = y_range[::-1] / nautical_miles_conversion  # Convert to nautical miles
#     percent_diff_at_x_center = percent_diff_at_x_center[::-1]

#     # Find the indices for the specified nautical miles and store results
#     results = {'File': file_path}
#     for nm in nautical_miles:
#         idx = np.abs(y_range - nm).argmin()
#         results[f'{nm} nautical miles'] = percent_diff_at_x_center[idx]
    
#     results_list.append(results)

# # Create a DataFrame and save to CSV
# results_df = pd.DataFrame(results_list)
# output_file = 'percent_diff.csv'
# results_df.to_csv(output_file, index=False)