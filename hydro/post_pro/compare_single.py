import numpy as np
import os
import csv
import matplotlib.pyplot as plt

def read_csv(file_name, header_adjustment):
    w_val = None
    Kt_H = []
    Kr_H = []
    power = []
    
    script_dir = os.path.dirname(__file__)
    data_folder = os.path.join(script_dir, '..','data','data_indiv')
    file_path = os.path.join(data_folder, file_name)
    
    with open(file_path, mode='r') as file:
        reader = csv.reader(file)
        header = next(reader)
        
        num_bodies = (len(header) - header_adjustment) // 2
        Kt_H = [None] * num_bodies
        Kr_H = [None] * num_bodies
        power = [None] * num_bodies
        
        for i, row in enumerate(reader):
            if i == 0:  # Row index 0 (1st row)
                w_val = float(row[0])
                for j in range(num_bodies):
                    Kt_H[j] = float(row[1 + j])
                    Kr_H[j] = float(row[1 + num_bodies + j])
                    power[j] = float(row[1 + 2*num_bodies + j])
                break
                
    return w_val, Kt_H, Kr_H, power

def process_files(file_names, header_adjustments):
    Kt_H_values = []
    Kr_H_values = []
    for file_name, header_adjustment in zip(file_names, header_adjustments):
        _, Kt_H, Kr_H, _ = read_csv(file_name, header_adjustment)
        Kt_H_values.append(np.array(Kt_H))
        Kr_H_values.append(np.array(Kr_H))
    return Kt_H_values, Kr_H_values

def calculate_percent_diff(values1, values2):
    return 100 * ((values2 - values1) / values2)

def plot_bars(data, x_positions, width, colors):
    bars = []
    for i, (values, color) in enumerate(zip(data, colors)):
        bars.append(plt.bar(x_positions + i * width, values, width, color=color))
    return bars

def annotate_bars(bars, percent_diffs):
    for i, (bar_group, pd_group) in enumerate(zip(bars, percent_diffs)):
        for j, bar in enumerate(bar_group):
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width() / 2.0, height, f'{pd_group[j - 1]:.2f}%', ha='center', va='bottom', fontsize=14, fontweight='bold')

# File names and header adjustments
file_names1 = ['PA_un_125.csv', 'PA_damp_125.csv', 'PA_react_125.csv']
file_names2 = ['OS_un_125.csv', 'OS_damp_125.csv', 'OS_react_125.csv']
file_names3 = ['atten_un_125.csv', 'atten_damp_125.csv', 'atten_react_125.csv']
file_names4 = ['break_125.csv']
header_adjustments = [1, 1, 1, 1]  # Adjust to include the fourth file

# Process files
Kt_H_values1, Kr_H_values1 = process_files(file_names1, header_adjustments[:3])
Kt_H_values2, Kr_H_values2 = process_files(file_names2, header_adjustments[:3])
Kt_H_values3, Kr_H_values3 = process_files(file_names3, header_adjustments[:3])
Kt_H_values4, Kr_H_values4 = process_files(file_names4, header_adjustments[3:])

# Calculate percentage differences for Kt_H and Kr_H
percent_diff_kt1 = [calculate_percent_diff(Kt_H_values1[0], Kt_H_values1[i]) for i in range(1, len(Kt_H_values1))]
percent_diff_kt2 = [calculate_percent_diff(Kt_H_values2[0], Kt_H_values2[i]) for i in range(1, len(Kt_H_values2))]
percent_diff_kt3 = [calculate_percent_diff(Kt_H_values3[0], Kt_H_values3[i]) for i in range(1, len(Kt_H_values3))]
# No percentage difference for the single file in the fourth group
percent_diff_kr1 = [calculate_percent_diff(Kr_H_values1[0], Kr_H_values1[i]) for i in range(1, len(Kr_H_values1))]
percent_diff_kr2 = [calculate_percent_diff(Kr_H_values2[0], Kr_H_values2[i]) for i in range(1, len(Kr_H_values2))]
percent_diff_kr3 = [calculate_percent_diff(Kr_H_values3[0], Kr_H_values3[i]) for i in range(1, len(Kr_H_values3))]
# No percentage difference for the single file in the fourth group

# Plotting the Kt_H values as a bar plot
num_bodies = len(Kt_H_values1[0])
x_labels_kt = [f'PA' for i in range(num_bodies)]
x_labels_kr = [f'OS' for i in range(num_bodies)]
x_labels_atten = [f'Atten' for i in range(num_bodies)]
x_labels_break = ['Break']  # Single bar label for the fourth file
x_labels = x_labels_kt + x_labels_kr + x_labels_atten + x_labels_break

# Create bar plots
width = 0.275  # Width of the bars
x = np.arange(num_bodies * 3 + 1)  # Adjust x to accommodate all bars, including the single bar

# Kt_H Plot
plt.figure(figsize=(20, 12))
colors = ['#E69F00', '#56B4E9', '#009E73']

bars_kt1 = plot_bars(Kt_H_values1, x[:num_bodies], width, colors)
bars_kt2 = plot_bars(Kt_H_values2, x[num_bodies:num_bodies*2], width, colors)
bars_kt3 = plot_bars(Kt_H_values3, x[num_bodies*2:num_bodies*3], width, colors)
bars_kt4 = plot_bars([Kt_H_values4[0]], x[num_bodies*3:], width, ['#E69F00'])

# Annotate bars with percentage differences for Kt_H
#annotate_bars(bars_kt1[1:], percent_diff_kt1)
#annotate_bars(bars_kt2[1:], percent_diff_kt2)
#annotate_bars(bars_kt3[1:], percent_diff_kt3)

# Get the height of the "break_125" bar
break_kt_height = Kt_H_values4[0][0]

plt.axhline(y=break_kt_height, color='gray', linestyle='--', linewidth=4,label='_nolegend_')

plt.ylim(0, 1.0)
plt.xticks(x, x_labels, fontsize=40, fontweight='bold', rotation=45)
plt.yticks(fontsize=30, fontweight='bold')
plt.ylabel('$K_t$', fontsize=45, fontweight='bold')
plt.legend(['Uncontrolled', 'Damped', 'Reactive'], loc='upper left', bbox_to_anchor=(1, 1),fontsize=40)
plt.tight_layout()
plt.savefig('kt_comparison.pdf')
plt.close()

# Kr_H Plot
plt.figure(figsize=(16, 10))
colors = ['#E69F00', '#56B4E9', '#009E73']

# Plotting Kr_H values for the first group
bars_kr1 = plot_bars(Kr_H_values1, x[:num_bodies], width, colors)

# Plotting Kr_H values for the second group
bars_kr2 = plot_bars(Kr_H_values2, x[num_bodies:num_bodies*2], width, colors)

# Plotting Kr_H values for the third group
bars_kr3 = plot_bars(Kr_H_values3, x[num_bodies*2:num_bodies*3], width, colors)

# Plotting Kr_H values for the fourth group
bars_kr4 = plot_bars([Kr_H_values4[0]], x[num_bodies*3:], width, ['#E69F00'])

# Annotate bars with percentage differences for Kr_H
annotate_bars(bars_kr1[1:], percent_diff_kr1)
annotate_bars(bars_kr2[1:], percent_diff_kr2)
annotate_bars(bars_kr3[1:], percent_diff_kr3)

# Get the height of the "break_125" bar for Kr_H
break_kr_height = Kr_H_values4[0][0]

plt.axhline(y=break_kr_height, color='gray', linestyle='--', linewidth=4,label='_nolegend_')

plt.ylim(0, 1.0)
plt.xticks(x, x_labels, fontsize=35, fontweight='bold', rotation=45)
plt.yticks(fontsize=30, fontweight='bold')
plt.ylabel('$K_r$', fontsize=45, fontweight='bold')
plt.legend(['Uncontrolled', 'Damped', 'Reactive'], loc='upper right', fontsize=25)
plt.tight_layout()
plt.savefig('kr_comparison.pdf')
plt.close()

print('kt_comparison.pdf saved')
print('kr_comparison.pdf saved')