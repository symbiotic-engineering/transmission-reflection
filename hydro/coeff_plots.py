import os
import csv
import matplotlib.pyplot as plt

def read_csv(file_name):
    w_vals = []
    Kt_H = []
    Kr_H = []
    
    script_dir = os.path.dirname(__file__)
    data_folder = os.path.join(script_dir, 'data')
    file_path = os.path.join(data_folder, file_name)
    
    with open(file_path, mode='r') as file:
        reader = csv.reader(file)
        header = next(reader)
        
        num_bodies = (len(header) - 1) // 2
        Kt_H = [[] for _ in range(num_bodies)]
        Kr_H = [[] for _ in range(num_bodies)]
        
        for row in reader:
            w_vals.append(float(row[0]))
            for i in range(num_bodies):
                Kt_H[i].append(float(row[1 + i]))
                Kr_H[i].append(float(row[1 + num_bodies + i]))
                
    return w_vals, Kt_H, Kr_H

def plot_data(w_vals, Kt_H, Kr_H, cud_colors, linestyles, label_prefix, marker_offset, color_offset, linestyle_offset):
    markers_list = ['o', 's', 'd', '^', 'p', '*', '<', 'p']  # List of markers for each line

    # Plot Kt lines
    for i, kt_h_values in enumerate(Kt_H):
        plt.plot(w_vals, kt_h_values, marker=markers_list[(i + marker_offset) % len(markers_list)], markersize=10, 
                 label=f'{label_prefix} $K_t$',
                 color=cud_colors[(i + color_offset) % len(cud_colors)], 
                 linestyle=linestyles[(i + linestyle_offset) % len(linestyles)], linewidth=3)
    
    # Plot Kr lines
    for i, kr_h_values in enumerate(Kr_H):
        plt.plot(w_vals, kr_h_values, marker=markers_list[(i + 1 + marker_offset) % len(markers_list)], markersize=10, 
                 label=f'{label_prefix} $K_r$',
                 color=cud_colors[(i + 1 + color_offset) % len(cud_colors)], 
                 linestyle=linestyles[(i + 1 + linestyle_offset) % len(linestyles)], linewidth=3)

# File names (adjust these according to your actual file names)
file_name1 = 'break_uncont.csv'
#file_name3 = 'atten_reactive.csv'
#file_name2 = 'atten_damp.csv'

# Read data from CSV files
w_vals1, Kt_H1, Kr_H1 = read_csv(file_name1)
#w_vals2, Kt_H2, Kr_H2 = read_csv(file_name2)
#w_vals3, Kt_H3, Kr_H3 = read_csv(file_name3)

# Plotting parameters
cud_colors = ['#E69F00', '#56B4E9', '#009E73', '#D55E00', '#CC79A7', '#000000']
linestyles = ['-', '--', ':', '-.', '-', '--', ':', '-.']

# Plot data with different labels and ensure unique colors and markers for each dataset
plot_data(w_vals1, Kt_H1, Kr_H1, cud_colors, linestyles, 'Uncontrolled', marker_offset=0, color_offset=0, linestyle_offset=0)
#plot_data(w_vals2, Kt_H2, Kr_H2, cud_colors, linestyles, 'Damped', marker_offset=2, color_offset=2, linestyle_offset=2)
#plot_data(w_vals3, Kt_H3, Kr_H3, cud_colors, linestyles, 'Reactive', marker_offset=4, color_offset=4, linestyle_offset=4)

# plt.legend(fontsize=15, markerscale=1)#, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('$\omega$ [rad/s]', fontsize=20)
plt.ylabel('Coefficient Value', fontsize=20)
plt.tight_layout()

print('tee')
plt.savefig('break_single.pdf')
print('hee')

