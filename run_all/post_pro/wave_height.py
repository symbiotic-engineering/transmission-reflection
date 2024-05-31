import numpy as np
import matplotlib.pyplot as plt

# Constants
xgrid = 3000
ygrid = 5000
mxc = 300
myc = 300
x_conversion = xgrid / (mxc + 1)
y_conversion = ygrid / (myc + 1)
x_investigated = int(1550 / x_conversion)       # dead x-center of array
y_investigated = int(4385 / y_conversion)       # y = 0 m from array

# Load and reshape data
def load_and_reshape(file_path):
    data = np.loadtxt(file_path, delimiter=',')
    return data.reshape(mxc + 1, myc + 1)

file_paths = [
    'blank_elevation.csv',
    'atten_1a.csv','atten_1b.csv','atten_1c.csv','atten_1d.csv','atten_2a.csv','atten_3b.csv',
    'break_1a.csv','break_1b.csv','break_1c.csv','break_1d.csv','break_2a.csv','break_3b.csv',
    'OS_1a.csv','OS_1b.csv','OS_1c.csv','OS_1d.csv','OS_2a.csv','OS_3b.csv',
    'PA_1a.csv','PA_1b.csv','PA_1c.csv','PA_1d.csv','PA_2a.csv','PA_3b.csv'
]

data_arrays = [load_and_reshape(f'/mnt/c/Users/ov162/transmission-reflection/data/{fp}') for fp in file_paths]

# Extract midline heights
midline_heights = [data[:y_investigated, x_investigated] for data in data_arrays]

# Calculate percent decrease
midline_height_blank = midline_heights[0]
percent_decrease = 100 * (midline_height_blank - np.array(midline_heights[1:])) / midline_height_blank

# Calculate distance from array
distance_from_array = abs(np.linspace(0, y_investigated * y_conversion, len(midline_heights[1])))# - ygrid)
#print(distance_from_array)
# Find the indices for the specified x-distances
x_distance_1nm = 5000 - 1852
x_distance_half_nm = 5000-(3704)
x_index_1nm = np.argmin(np.abs(distance_from_array - x_distance_1nm))
x_index_half_nm = np.argmin(np.abs(distance_from_array - x_distance_half_nm))

# Extract values at the specified x-distances for each file
percent_decrease_1nm_values = percent_decrease[:, x_index_1nm]
percent_decrease_half_nm_values = percent_decrease[:, x_index_half_nm]

# Print the extracted values for each file
for i, file_path in enumerate(file_paths[1:]):  # Skipping the blank file
    #print(f"{file_path} - Percent decrease at 1852 meters: {percent_decrease_1nm_values[i]:.2f}%")
    print(f"{file_path} - Percent decrease at 926 meters: {percent_decrease_half_nm_values[i]:.2f}%")


# Define color-blind friendly colors
cud_colors = ['#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7', '#000000', '#8B4513']

# Plot function
def plot_data(x, y_data, ylabel, filename):
    plt.figure()
    for y, label, color, linestyle in zip(y_data, plot_labels, cud_colors, ['-', '--', '-.', ':','-','--','-.',':']):
        plt.plot(x, y, label=label, color=color, linestyle=linestyle)
    plt.axvline(x=1852, color='black', linestyle=':', label='1 nm')
    plt.axvline(x=(1852)/2, color='black', linestyle='-')
    plt.xlabel('Distance from Array [m]')
    plt.ylabel(ylabel)
    plt.legend()
    #plt.savefig(f'/mnt/c/Users/ov162/transmission-reflection/run_all/post_pro/{filename}.pdf')
    plt.show()

# Plot wave heights
plot_labels = ['Breakwater', 'PA', 'Attenuator','OS']
#plot_data(distance_from_array, midline_heights[1:], 'Wave Height [m]', 'compare_bodies')

# Plot percent decrease
plot_data(distance_from_array, percent_decrease, 'Decrease in Wave Height [%]', '1d_percent_decrease')
