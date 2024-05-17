import numpy as np
import matplotlib.pyplot as plt

# Constants
xgrid = 3000
ygrid = 5000
mxc = 300
myc = 300
x_conversion = xgrid / (mxc + 1)
y_conversion = ygrid / (myc + 1)
x_positions = [1490, 1590, 1690, 1420, 1520, 1620]
ya = 4500
yb = 4400
H = 1.3832
x_investigated = int(1550 / x_conversion)
y_investigated = int(4385 / y_conversion)

# Load and reshape data
def load_and_reshape(file_path):
    data = np.loadtxt(file_path, delimiter=',')
    return data.reshape(mxc + 1, myc + 1)

file_paths = [
    'blank_elevation.csv',
    # 'PA_elevation.csv',
    # 'OSWEC_elevation.csv',
    # 'OSWEC_3bodcoeff.csv',
    # 'atten_elevation.csv',
    # 'PA_3bodcoeff.csv',
    # 'atten_3bodcoeff.csv',
    # 'break_elevation.csv',
    # 'break_3bodcoeff.csv',
    # 'PA6_elevation.csv',
    # 'PA_6bodcoeff.csv',
    'OSWEC6_elevation.csv',
    'OSWEC_6bodcoeff.csv',
    'atten6_elevation.csv',
    'atten_6bodcoeff.csv'
    # 'break6_elevation.csv',
    # 'break_6bodcoeff.csv'
]

data_arrays = [load_and_reshape(f'/mnt/c/Users/ov162/transmission-reflection/data/{fp}') for fp in file_paths]

# Extract midline heights
midline_heights = [data[:y_investigated, x_investigated] for data in data_arrays]

# Calculate percent decrease
midline_height_blank = midline_heights[0]
percent_decrease = 100 * (midline_height_blank - np.array(midline_heights[1:])) / midline_height_blank

# Calculate distance from array
distance_from_array = abs(np.linspace(0, y_investigated * y_conversion, len(midline_heights[1])) - ygrid)

# Define color-blind friendly colors
cud_colors = ['#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7', '#000000', '#8B4513']

# Plot function
def plot_data(x, y_data, ylabel, filename):
    plt.figure()
    for y, label, color, linestyle in zip(y_data, plot_labels, cud_colors, ['-', '--', '-.', ':','-','--','-.',':']):
        plt.plot(x, y, label=label, color=color, linestyle=linestyle)
    #plt.axvline(x=1852, color='black', linestyle='-', label='nm')
    #plt.axvline(x=1852*2, color='black', linestyle='-')
    plt.xlabel('Distance from Array [m]')
    plt.ylabel(ylabel)
    plt.legend()
    plt.savefig(f'/mnt/c/Users/ov162/transmission-reflection/run_all/post_pro/{filename}.pdf')
    plt.show()

# Plot wave heights
# plot_labels = ['PA', 'OSWEC', 'Attenuator', 'Breakwater', 'PA6', 'OSWEC6','atten6','break6']
plot_labels = ['OSWEC', 'OSWEC Indiv', 'Attenuator','Attenuator Indiv']
plot_data(distance_from_array, midline_heights[1:], 'Wave Height [m]', 'compare_bodies')

# Plot percent decrease
plot_data(distance_from_array, percent_decrease, 'Decrease in Wave Height [%]', 'percent_decrease_compare_bodies')
