import body
import solve
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

xtrans = np.array([0,0]) 
ytrans = np.array([50,50]) 
farm = False
controls = False
attenuator = False
B = 0
depth = 40
w = 1.25
rad = True
N = 3
D = 30

# for PA
r = 10              # radius [m]
l = 5               # length [m]
points = np.linspace(1,30,30)

# Initialize lists to store results
panels = []
tot_el_1pt = []

for i in range(len(points)):
# Compute cosine spacing for each dimension
    nr = int((r/l)*points[i])
    print(nr)
    nz = int((l/r)*points[i])
    ntheta = int((2*np.pi)*points[i])

    array, rel_dim = body.PA(xtrans, ytrans, farm,nr,ntheta,nz)
    #array, rel_dim = body.attenuator(xtrans, ytrans, farm, D,nr,ntheta,nz)
    diff_result, rad_result, RAO_vals, lam = solve.hydro(array, B, depth, w, farm, controls)
    
    # Calculate elevation (replace with actual implementation)
    # total, incoming_fse, x1, x2, nx, y1, y2, ny = solve.elevation(res, lam, diff_result, rad_result, RAO_vals, farm, rad, controls, N, attenuator, rel_dim)
    
    # Extract the specific point from the array 'total'
    specific_point = RAO_vals[0, 0] # np.mean(total)
    
    # Append the extracted value to tot_el_1pt list
    tot_el_1pt.append(specific_point)
    
    # Append the radial mesh value
    panels.append(points)
df = pd.DataFrame(tot_el_1pt, columns=['TotEl1pt'])
# Define the CSV file path
csv_file_path = 'PA_RAO.csv'
# Write the DataFrame to a CSV file
df.to_csv(csv_file_path, index=False)

# Calculate the percent difference between consecutive tot_el_1pt values
percent_diff = [0]  # Start with 0 for the first value as there's no previous value to compare
for i in range(1, len(tot_el_1pt)):
    diff = 100 * (tot_el_1pt[i] - tot_el_1pt[i - 1]) / tot_el_1pt[i - 1]
    percent_diff.append(diff)

# Plot the percent difference
plt.figure(figsize=(12, 6))
plt.plot(panels, percent_diff, marker='o',color='red')
plt.xlabel('Z-Direction Number of Panels',fontsize=20)
plt.ylabel('Percent Difference [%]',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend()
plt.tight_layout()
plt.savefig('PA_panel_conv.pdf')
plt.show()
