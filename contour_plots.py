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
    'PA_large.csv'
]

# Generate x and y coordinates
x = np.linspace(0, xgrid, mxc + 1)
y = np.linspace(0, ygrid, myc + 1)

# Create contour plots
for file_name in file_paths:
    # Load data
    file_path = os.path.join(data_folder, file_name)
    data = load_and_reshape(file_path)

    x1, x2, x3, x4 = 1490,1550,1460,1520
    x_pos = [x1,x2,x3,x4,
         x1 + 200, x2 + 200, x3 + 200, x4 + 200,
         x1 - 200, x2 - 200, x3 - 200, x4 - 200,
         x1 + 400, x2 + 400, x3 + 400, x4 + 400,
         x1 - 400, x2 - 400, x3 - 400, x4 - 400]
    ya, yb = 4500, 4460

    # Calculate wecx and wecy for the scatter plot
    wecx = [i + 10 for i in x_pos]
    wecy = [ya,ya,yb,yb,
            ya,ya,yb,yb,
            ya,ya,yb,yb,
            ya,ya,yb,yb,
            ya,ya,yb,yb]
    
    # Plot contour and scatter
    plt.contourf(x, y, data, levels=50, cmap='viridis')
    plt.scatter(wecx, wecy, marker='o', color='red', s=3, linewidth=2)
    cbar = plt.colorbar()                                                   # Create the colorbar
    cbar.set_label('Wave Height [m]', fontsize=16)                          # Set the label with the desired font size
    cbar.ax.tick_params(labelsize=14)                                       # Set the font size for the colorbar ticks
    plt.xlabel('x [m]', fontsize=20)
    #plt.ylabel('y [m]', fontsize=20)
    plt.xticks(fontsize=15)
    #plt.yticks(fontsize=15)
    ax = plt.gca()
    ax.tick_params(left=False, bottom=True, labelleft=False, labelbottom=True)
    #plt.text(1.025, 1.10, 'e', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top', ha='left')
    
    # Save plot as PDF
    plt.tight_layout()
    plt.savefig(f'{file_name[:-4]}_contour.pdf', format='pdf')
    plt.close()

print("Contour plots created and saved as PDF files.")