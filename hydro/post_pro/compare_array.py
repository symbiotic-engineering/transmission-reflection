import numpy as np
import os
import csv
import matplotlib.pyplot as plt

def read_csv(file_name):
    Kt_H = []
    Kr_H = []
    
    # Navigate up two levels from the post_pro folder to transmission-reflection, then into the data folder
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, '..'))
    data_folder = os.path.join(project_root, 'data')
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
    
    # Navigate up two levels from the post_pro folder to transmission-reflection, then into the data folder
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, '..'))
    data_folder = os.path.join(project_root, 'data')
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
oswec_array = ['OS_stag_uncont.csv']
pa_array = ['PA_stag_uncont.csv']
atten_array = ['atten_stag_uncont.csv']
break_array = ['break_stag.csv']
oswec_iso = 'OS_uncont.csv'
pa_iso = 'PA_uncont.csv'
atten_iso = 'atten_uncont.csv'
break_iso = 'break.csv'

# Process files
Kt_OS_array, Kr_OS_array = process_files(oswec_array)
Kt_PA_array, Kr_PA_array = process_files(pa_array)
Kt_atten_array, Kr_atten_array = process_files(atten_array)
Kt_break_array, Kr_break_array = process_files(break_array)

Kt_OS, Kr_OS = read_single_file(oswec_iso)
Kt_PA, Kr_PA = read_single_file(pa_iso)
Kt_atten, Kr_atten = read_single_file(atten_iso)
Kt_break, Kr_break = read_single_file(break_iso)

# Define your arrays and corresponding file names in a dictionary
data_dict = {
    'OS': {'array': oswec_array, 'iso_file': oswec_iso},
    'PA': {'array': pa_array, 'iso_file': pa_iso},
    'attenuator': {'array': atten_array, 'iso_file': atten_iso},
    'breakwater': {'array': break_array, 'iso_file': break_iso}
}

# Initialize dictionaries to store the results
Kt_array_dict = {}
Kr_array_dict = {}
Kt_iso_dict = {}
Kr_iso_dict = {}
percent_difference_dict = {}

# Loop through the dictionary
for key, value in data_dict.items():
    # Process the files for each case
    Kt_array_dict[key], Kr_array_dict[key] = process_files(value['array'])
    
    # Read the corresponding ISO file
    Kt_iso_dict[key], Kr_iso_dict[key] = read_single_file(value['iso_file'])
    
    # Calculate the percent difference
    percent_difference_dict[key] = [
        abs(100 * (Kt_iso_dict[key] - coeff) / Kt_iso_dict[key]) for coeff in Kt_array_dict[key]
    ]
    
    # Print the percent difference for each case
    print(f"{key} percent difference: {percent_difference_dict[key]}")

# Colors and markers for different cases
colors = {'OS': '#E69F00', 'PA': '#56B4E9', 'attenuator': '#009E73', 'breakwater': '#D55E00'}
markers = {'OS': 'o', 'PA': 's', 'attenuator': '^', 'breakwater': '*'}

# Labels for each body
bodies = ['Body 1', 'Body 2', 'Body 3']

plt.figure(figsize=(7, 6))

# Loop through each body and plot the percent differences
for i in range(3):  # Loop over the three bodies
    for key in percent_difference_dict:  # Loop over OS, PA, and attenuator
        plt.scatter(i + 1, percent_difference_dict[key][0][i], 
                    color=colors[key], marker=markers[key], s=250,
                    label=f'{key}' if i == 0 else "")

# Customize the plot
#plt.xticks([1, 2, 3], bodies, fontsize=20)
plt.yticks(fontsize=20)
#plt.ylabel('$K_t$ Change [%]',fontsize=20)
plt.ylim([0,80])
ax = plt.gca()
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)
#ax.spines['left'].set_visible(False)
#ax.spines['bottom'].set_visible(False)

ax.xaxis.set_major_locator(plt.MultipleLocator(0.25))

ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
#plt.legend(loc='center',fontsize=20)
plt.grid(True)
plt.text(1.025, 1.05, 'b', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top', ha='left')

# Show and save the plot
plt.savefig('percent_difference_plot.pdf')
plt.show()


# # Plotting the Kt_H values as a bar plot
# x_labels = ['Body 1', 'Body 2', 'Body 3']

# # Create bar plots
# width = 0.2  # Width of the bars
# x = np.arange(len(x_labels))

# # Offsets
# offset = width/3  # Offset for additional file's data
# colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442']

# # Kt_H Plot
# #plt.figure(figsize=(18, 12))
# plt.figure(figsize=(14, 12))
# plt.bar(x - width / 2, Kt_H_values1[0], width, color=colors[1],label='Uncont')
# plt.bar(x + width / 2, Kt_H_values2[0], width, color=colors[2],label='Damp')
# plt.bar(x + width*1.5, Kt_H_values3[0], width, color=colors[3],label='Reactive')
# plt.bar(x - width / 2, Kt_H_single[0], width, facecolor='none', linewidth=4,hatch='/', label='Isolated', edgecolor='black')
# plt.bar(x + width / 2, Kt_H_single2[0], width, facecolor='none', linewidth=4,hatch='/', edgecolor='black')
# plt.bar(x + width*1.5, Kt_H_single3[0], width, facecolor='none', linewidth=4,hatch='/', edgecolor='black')

# plt.ylim(0, 1.2)
# plt.axhline(y=1.0, color='gray', linestyle='--', linewidth=4,label='_nolegend_')
# plt.xticks(x, x_labels, fontsize=40, fontweight='bold', rotation=45)
# plt.yticks(fontsize=40, fontweight='bold')
# plt.ylabel('$K_t$', fontsize=50, fontweight='bold')
# #plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=40)
# plt.tight_layout()
# plt.savefig('kt_comp_reg.pdf')
# plt.close()

# # Kr_H Plot
# plt.figure(figsize=(18, 12))
# plt.bar(x - width / 2, Kr_H_values1[0], width, color=colors[1], label='Uncont')
# plt.bar(x + width / 2, Kr_H_values2[0], width, color=colors[2], label='Damp')
# plt.bar(x + width * 1.5, Kr_H_values3[0], width, color=colors[3], label='Reactive')
# plt.bar(x - width / 2, Kr_H_single[0], width, facecolor='none', linewidth=4, hatch='/', label='Isolated', edgecolor='black')
# plt.bar(x + width / 2, Kr_H_single2[0], width, facecolor='none', linewidth=4, hatch='/', edgecolor='black')
# plt.bar(x + width * 1.5, Kr_H_single3[0], width, facecolor='none', linewidth=4, hatch='/', edgecolor='black')

# # Add similar styling
# plt.ylim(0, 1.2)
# plt.axhline(y=1.0, color='gray', linestyle='--', linewidth=4, label='_nolegend_')
# plt.xticks(x, x_labels, fontsize=35, fontweight='bold', rotation=45)
# plt.yticks(fontsize=35, fontweight='bold')
# plt.ylabel('$K_r$', fontsize=45, fontweight='bold')
# plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=35)
# plt.tight_layout()
# plt.savefig('kr_comp_reg.pdf')
# plt.close()

# print('kt_comp_array.pdf saved')
# print('kr_comp_array.pdf saved')