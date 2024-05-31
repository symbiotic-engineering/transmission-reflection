import body
import solve
import numpy as np
import matplotlib.pyplot as plt

xtrans = np.array([0,0]) 
ytrans = np.array([50,50]) 
farm = False
controls = False
attenuator = False
B = 0
depth = 40
w = 0.7
rad = True
N = 3
D = 30

# for the point absorber
# nz_values = np.linspace(1, 50, 49)
# ntheta_values = np.linspace(5,30,24)
# nr_values = np.linspace(5,50,44)

# for OSWEC
# nt = 3         # this might need to be 8... eek!
# nh = 10        # this was good!        
# nw = 10         # this was good!
# nw_values = np.linspace(1,30,29)#10         # number of panels along width (x)
# nt_values = np.linspace(1,20,18)          # number of panels along thickness (y)
# nh_values = np.linspace(1,30,29)         # number of panels along height (z)

# for breakwater
# nw = 10         # number of panels along width (x)
# nt = 5          # number of panels along thickness (y) might need to be 13, eek!
# nh = 2          # number of panels along height (z) should probably be 20 yikes!!
# nw_values = np.linspace(1,30,29)         # number of panels along width (x)
# nt_values = np.linspace(1,20,18)          # number of panels along thickness (y)
# nh_values = np.linspace(1,30,29)         # number of panels along height (z)

# for attenuator
# nr = 2        # 10 or 15 would be better
# ntheta = 5     # 12 would be better
# nz = 10            # 10 was perfect!
# nz_values = np.linspace(1, 50, 49)
# ntheta_values = np.linspace(3,40,39)
# nr_values = np.linspace(1,20,19)

# Initialize lists to store results
radial_mesh = []
tot_el_1pt = []
res_values = np.linspace(1,4,4)

for i in range(len(res_values)):
    res = int(res_values[i])
    array, rel_dim = body.PA(xtrans, ytrans, farm)
    #array, rel_dim = body.attenuator(xtrans, ytrans, farm, D,nr,ntheta,nz)
    diff_result, rad_result, RAO_vals, lam = solve.hydro(array, B, depth, w, farm, controls)
    
    # Calculate elevation (replace with actual implementation)
    total, incoming_fse, x1, x2, nx, y1, y2, ny = solve.elevation(res, lam, diff_result, rad_result, RAO_vals, farm, rad, controls, N, attenuator, rel_dim)
    
    # Extract the specific point from the array 'total'
    specific_point = np.mean(total)# RAO[0, 0]
    
    # Append the extracted value to tot_el_1pt list
    tot_el_1pt.append(specific_point)
    
    # Append the radial mesh value
    radial_mesh.append(res)

# Calculate the percent difference between consecutive tot_el_1pt values
percent_diff = [0]  # Start with 0 for the first value as there's no previous value to compare
for i in range(1, len(tot_el_1pt)):
    diff = 100 * (tot_el_1pt[i] - tot_el_1pt[i - 1]) / tot_el_1pt[i - 1]
    percent_diff.append(diff)

# Plot the percent difference
plt.figure(figsize=(12, 6))
plt.plot(radial_mesh, percent_diff, marker='o',color='red')
plt.xlabel('Resolution Multiplier to Wavelength',fontsize=17)
plt.ylabel('Percent Difference [%]',fontsize=17)
#plt.title('Percent Difference in Total Elevation vs Z-direction Mesh Values')
plt.legend()
plt.tight_layout()
plt.savefig('07_res_smaller_conv_perc.pdf')
plt.show()
