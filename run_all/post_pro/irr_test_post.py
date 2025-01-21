import numpy as np
import matplotlib.pyplot as plt

# Grid and physical dimensions
xgrid = 3000  # Physical grid width in meters
ygrid = 5000  # Physical grid height in meters
mxc = 300     # Number of grid cells in x
myc = 500     # Number of grid cells in y

# Conversion factors from physical coordinates to indices
x_conversion = xgrid / (mxc + 1)
y_conversion = ygrid / (myc + 1)

# Convert physical coordinates to grid indices
x_investigated = int(1550 / x_conversion)  # Index for x = 1550
y_end = int(4482 / y_conversion)           # Index for y = 4482

# Specific y-coordinates to investigate
y_specific_1 = 4482 - 1852
y_specific_2 = 4482 - (1852 / 2)

y_specific_index_1 = int(y_specific_1 / y_conversion)
y_specific_index_2 = int(y_specific_2 / y_conversion)

# Load and reshape data
def load_and_reshape(file_path):
    data = np.loadtxt(file_path, delimiter=',')
    return data.reshape(myc + 1, mxc + 1)  # Reshape to match (myc + 1, mxc + 1)

# File paths
blank_path = 'blank.csv'
irregular_path = 'irregular_test.csv'
regular_path = 'OSr_1a.csv'

# Load data
blank_data = load_and_reshape(blank_path)
irregular_data = load_and_reshape(irregular_path)
regular_data = load_and_reshape(regular_path)

# Extract centerline data
centerline_blank = blank_data[:y_end + 1, x_investigated]  # Data from y=0 to y=4482
centerline_irregular = irregular_data[:y_end + 1, x_investigated]
centerline_regular = regular_data[:y_end + 1, x_investigated]

# Calculate percent difference
percent_difference = 100 * (centerline_blank - centerline_irregular) / centerline_blank
reg_percent_difference = 100 * (centerline_blank - centerline_regular) / centerline_blank

# Generate y-axis (distance) values
y_values = np.linspace(0, y_end * y_conversion, len(centerline_blank))

# Extract and display values at specific y-coordinates
value_blank_1 = centerline_blank[y_specific_index_1]
value_irregular_1 = centerline_irregular[y_specific_index_1]
percent_difference_1 = percent_difference[y_specific_index_1]

value_regular_1 = centerline_regular[y_specific_index_1]
reg_percent_difference_1 = reg_percent_difference[y_specific_index_1]

value_blank_2 = centerline_blank[y_specific_index_2]
value_irregular_2 = centerline_irregular[y_specific_index_2]
percent_difference_2 = percent_difference[y_specific_index_2]

value_regular_2 = centerline_regular[y_specific_index_2]
reg_percent_difference_2 = reg_percent_difference[y_specific_index_2]

print(f"At y = {y_specific_1} meters:")
print(f"  Percent difference: {percent_difference_1:.2f}%")

print(f"\nAt y = {y_specific_2} meters:")
print(f"  Percent difference: {percent_difference_2:.2f}%")

print(f"At y = {y_specific_1} meters:")
print(f"  Reg Percent difference: {reg_percent_difference_1:.2f}%")

print(f"\nAt y = {y_specific_2} meters:")
print(f"  Reg Percent difference: {reg_percent_difference_2:.2f}%")

plt.figure(figsize=(12, 8))  # Larger figure size for better visibility

# Define colorblind-friendly colors
color_irr = '#0072B2'  # Blue
color_reg = '#D55E00'  # Orange
color_1nm = '#CC79A7'  # Purple
color_0_5nm = '#009E73'  # Green

plt.plot(y_values / 1852, percent_difference, label='Irregular', color=color_irr, linewidth=10)
plt.plot(y_values / 1852, reg_percent_difference, label='Regular', color=color_reg, linestyle=':', linewidth=14)
plt.axvline(y_specific_1 / 1852, color=color_1nm, linestyle='--', linewidth=4, label='y = 1 nm')
plt.axvline(y_specific_2 / 1852, color=color_0_5nm, linestyle='--', linewidth=4, label='y = 0.5 nm')
plt.xlabel('Distance along y-axis [nm]', fontsize=24)
plt.ylabel('Wave Height Reduction [%]', fontsize=24)
plt.grid(True, linewidth=1.5, linestyle=':')
plt.legend(fontsize=22)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.tight_layout()
plt.savefig('compare_reg_irr.pdf', dpi=300)  # High resolution for better quality
plt.show()

