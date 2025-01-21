import numpy as np
import os
import csv
import matplotlib.pyplot as plt

def read_csv(file_name):
    Kt_H = []
    Kr_H = []
    
    script_dir = os.path.dirname(__file__)
    data_folder = os.path.join(script_dir, '..', 'data')
    file_path = os.path.join(data_folder, file_name)
    
    with open(file_path, mode='r') as file:
        reader = csv.reader(file)
        header = next(reader)
        
        for i, row in enumerate(reader):
            if i == 5:  # Index for w=1.25rad/s (sixth row, index 5)
                Kt_H = [float(row[1]), float(row[2]), float(row[3])]
                Kr_H = [float(row[4]), float(row[5]), float(row[6])]
                break
                
    return Kt_H, Kr_H

def read_single_file(file_name):
    Kt_H = []
    Kr_H = []
    
    script_dir = os.path.dirname(__file__)
    data_folder = os.path.join(script_dir, '..', 'data')
    file_path = os.path.join(data_folder, file_name)
    
    with open(file_path, mode='r') as file:
        reader = csv.reader(file)
        header = next(reader)
        
        for i, row in enumerate(reader):
            if i == 5:  # Index for w=1.25rad/s (sixth row, index 5)
                Kt_H = [float(row[1])]
                Kr_H = [float(row[2])]
                break
                
    return Kt_H, Kr_H

def process_files(file_names):
    Kt_H_values = []
    Kr_H_values = []
    for file_name in file_names:
        Kt_H, Kr_H = read_csv(file_name)
        Kt_H_values.append(np.array(Kt_H))
        Kr_H_values.append(np.array(Kr_H))
    return Kt_H_values, Kr_H_values

# File names
regular = ['OS_reg_uncont.csv']
staggered = ['OS_reg_damp.csv']
staggered3 = ['OS_reg_reactive.csv']
isolated1 = 'OS_uncont.csv'
isolated2 = 'OS_damp.csv'
isolated3 = 'OS_reactive.csv'

# Process files
Kt_H_values1, Kr_H_values1 = process_files(regular)
Kt_H_values2, Kr_H_values2 = process_files(staggered)
Kt_H_values3, Kr_H_values3 = process_files(staggered3)
Kt_H_single, Kr_H_single = read_single_file(isolated1)
Kt_H_single2, Kr_H_single2 = read_single_file(isolated2)
Kt_H_single3, Kr_H_single3 = read_single_file(isolated3)

# Plotting the Kt_H values as a bar plot
x_labels = ['Body 1', 'Body 2', 'Body 3']

# Create bar plots
width = 0.2  # Width of the bars
x = np.arange(len(x_labels))

# Offsets
offset = width/3  # Offset for additional file's data
colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442']

# Kt_H Plot
plt.figure(figsize=(18, 12))
plt.bar(x - width / 2, Kt_H_values1[0], width, color=colors[1],label='Uncont')
plt.bar(x + width / 2, Kt_H_values2[0], width, color=colors[2],label='Damp')
plt.bar(x + width*1.5, Kt_H_values3[0], width, color=colors[3],label='Reactive')
plt.bar(x - width / 2, Kt_H_single[0], width, facecolor='none', linewidth=4,hatch='/', label='Isolated', edgecolor='black')
plt.bar(x + width / 2, Kt_H_single2[0], width, facecolor='none', linewidth=4,hatch='/', edgecolor='black')
plt.bar(x + width*1.5, Kt_H_single3[0], width, facecolor='none', linewidth=4,hatch='/', edgecolor='black')
#plt.bar(x[1] + offset, Kt_H_single[0], width, facecolor='none', linewidth=4,hatch='/', label='Isolated', edgecolor='black')

plt.ylim(0, 1.2)
plt.axhline(y=1.0, color='gray', linestyle='--', linewidth=4,label='_nolegend_')
plt.xticks(x, x_labels, fontsize=35, fontweight='bold', rotation=45)
plt.yticks(fontsize=35, fontweight='bold')
plt.ylabel('$K_t$', fontsize=45, fontweight='bold')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=35)
plt.tight_layout()
plt.savefig('kt_comp_reg.pdf')
plt.close()

# Kr_H Plot
plt.figure(figsize=(18, 12))
plt.bar(x - width / 2, Kr_H_values1[0], width, color=colors[1],label='Regular')
plt.bar(x + width / 2, Kr_H_values2[0], width, color=colors[2],label='Staggered')
plt.bar(x[1] + offset, Kr_H_single[0], width, color=colors[3], hatch='/', label='Isolated', edgecolor='black')

plt.ylim(0, 1.0)
plt.xticks(x, x_labels, fontsize=35, fontweight='bold', rotation=45)
plt.yticks(fontsize=35, fontweight='bold')
plt.ylabel('$K_r$', fontsize=45, fontweight='bold')
plt.legend(loc='upper right', fontsize=35)
plt.tight_layout()
plt.savefig('kr_comp_reg.pdf')
plt.close()

print('kt_comp_array.pdf saved')
print('kr_comp_array.pdf saved')