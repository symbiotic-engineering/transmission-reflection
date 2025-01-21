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
data_folder = 'data/'

# Define file paths
file_paths = [
    #'blank.csv',
    #'OSr_1a.csv', 
    #'OSr_2a.csv',
    #'OSr_1b.csv', 
    #'OSr_3b.csv',
    #'OSr_1c.csv',
    #'OSr_1d.csv',
    #'attenr_1a.csv',
    #'attenr_2a.csv',
    #'attenr_1b.csv',
    #'attenr_3b.csv',
    #'attenr_1c.csv',
    #'attenr_1d.csv',
    #'break_1a.csv',
    #'break_2a.csv',
    #'break_1b.csv',
    #'break_3b.csv', 
    #'break_1c.csv',
    #'break_1d.csv',
    #'PAr_1a.csv',
    #'PAr_2a.csv',
    #'PAr_1b.csv',
    #'PAr_3b.csv',
    'PAr_1c.csv',
    #'PAr_1d.csv'
]

# Generate x and y coordinates
x = np.linspace(0, xgrid, mxc + 1)
y = np.linspace(0, ygrid, myc + 1)

# Create contour plots
for file_name in file_paths:
    # Load data
    file_path = os.path.join(data_folder, file_name)
    data = load_and_reshape(file_path)

    # Determine scatter plot coordinates based on the file name
    if file_name.endswith('a.csv'):
        x_pos = [1490, 1550, 1610] 
        ya, yb = 4500, 4500
    elif file_name.endswith('b.csv'):
        x_pos = [1490, 1550, 1610] 
        ya, yb = 4500, 4450 
    elif file_name.endswith('c.csv'):
        x_pos = [1290, 1350, 1410, 1470, 1530, 1590]
        ya, yb = 4500, 4500      
    elif file_name.endswith('d.csv'):
        x_pos = [1490, 1550, 1610, 1440, 1500, 1560]
        ya, yb = 4500, 4450

    # Calculate wecx and wecy for the scatter plot
    wecx = [i + 10 for i in x_pos]
    wecy = [yb if i > 2 else ya for i in range(len(wecx))]
    if file_name.endswith('b.csv'):
        wecy = [4450, 4500, 4450]
    
    # Plot contour and scatter
    plt.contourf(x, y, data, levels=50, cmap='viridis')
    plt.scatter(wecx, wecy, marker='_', color='red', s=25, linewidth=2)
    # cbar = plt.colorbar()  # Create the colorbar
    # cbar.set_label('Wave Height [m]', fontsize=16)  # Set the label with the desired font size
    # cbar.ax.tick_params(labelsize=14)  # Set the font size for the colorbar ticks
    plt.xlabel('X [m]', fontsize=20)
    plt.ylabel('Y [m]', fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    ax = plt.gca()
    ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    plt.text(1.025, 1.10, 'e', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top', ha='left')
    
    # Save plot as PDF
    plt.tight_layout()
    plt.savefig(f'{file_name[:-4]}_contour.pdf', format='pdf')
    plt.close()

print("Contour plots created and saved as PDF files.")