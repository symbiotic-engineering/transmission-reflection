import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
wind_dir = os.path.join(parent_dir, 'wind')
sys.path.append(wind_dir)
import morrison

# Constants
mxc = 300
myc = 700
xgrid = 3000  # Total width in meters
ygrid = 7000  # Total height in meters

# Load and reshape data
def load_and_reshape(file_path):
    data = np.loadtxt(file_path, delimiter=',')
    return data.reshape(myc + 1, mxc + 1)

# Set the correct data folder path
data_folder = '../data/'

# Define file paths
file_paths = [
    #'blank.csv',
    #'OSr_1a.csv', 
    'OS_SF_Maha.csv'
]

# Generate x and y coordinates
x = np.linspace(0, xgrid/1852, mxc + 1)
y = np.linspace(0, ygrid/1852, myc + 1)

# Create contour plots
for file_name in file_paths:
    # Load data
    file_path = os.path.join(data_folder, file_name)
    data = load_and_reshape(file_path)

    blank_filepath = os.path.join(data_folder, 'blank_Maha.csv')
    blank = load_and_reshape(blank_filepath)

    percent_data = (data/blank)
    percent_diff_D = morrison.fatigue_damage(percent_data)
    print('size of pdiff',np.size(percent_diff_D[0]))
    percent_diff_D.shape  = (myc + 1, mxc + 1)
    print('size of pdiff',np.size(percent_diff_D[0,:]))
    print('size of pdiff',np.size(percent_diff_D[0,:]))

    x1, x2, x3, x4 = 1490/1852,1550/1852,1460/1852,1520/1852
    x_pos = [x1,x2,x3,x4,
         x1 + 200/1852, x2 + 200/1852, x3 + 200/1852, x4 + 200/1852,
         x1 - 200/1852, x2 - 200/1852, x3 - 200/1852, x4 - 200/1852,
         x1 + 400/1852, x2 + 400/1852, x3 + 400/1852, x4 + 400/1852,
         x1 - 400/1852, x2 - 400/1852, x3 - 400/1852, x4 - 400/1852]
    ya, yb = 6500/1852, 6460/1852

    # Calculate wecx and wecy for the scatter plot
    wecx = [i + 2/1852 for i in x_pos]
    wecy = [ya,ya,yb,yb,
            ya,ya,yb,yb,
            ya,ya,yb,yb,
            ya,ya,yb,yb,
            ya,ya,yb,yb]

    # Plot contour and scatter
    plt.figure(figsize=(9.25,12))
    pcm = plt.pcolormesh(x, y, np.abs(percent_diff_D),vmin=0.00,vmax=70.00,cmap = 'viridis_r')    #,vmin=1.95, vmax=2.25)#, levels=100, cmap='viridis')
    pcm.set_edgecolor("face")
    cbar = plt.colorbar()                                                   # Create the colorbar
    cbar.set_label('Fatigue Damage Reduction [%]', fontsize=20,rotation=270,labelpad=40)                          # Set the label with the desired font size
    cbar.ax.tick_params(labelsize=18)                                       # Set the font size for the colorbar ticks
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.scatter(wecx, wecy, marker='_', color='red', s=10, linewidth=2)
    plt.xlabel('x [nm]', fontsize=20)
    #plt.ylabel('y [nm]', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    ax = plt.gca()
    ax.tick_params(left=False, bottom=True, labelleft=False, labelbottom=True)
    #plt.text(1.025, 1.10, 'e', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top', ha='left')
    
    # Save plot as PDF
    plt.tight_layout()
    plt.savefig(f'{file_name[:-4]}_testcontour.pdf', format='pdf')
    plt.close()

print("Contour plots created and saved as PDF files.")