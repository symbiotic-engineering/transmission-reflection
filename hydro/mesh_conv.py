import body
import solve
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

xtrans = [0,0]
ytrans = [-50,50]
farm = False
controls = False
attenuator = False
B = 0
depth = 40
w1 = 1.0
w2 = 1.25
w3 = 1.3
rad = True
N = 3
D = 30

# Initialize lists to store results
panels = []
RAO_07 = []
RAO_125 = []
RAO_13 = []
points = np.linspace(0.5,6.0,12)

for i in range(len(points)):
    res = points[i]
    print('res',res)
#     nt = int((th/(2*h))*points[i])
#     nh = int(points[i])
#     if nh < 1:
#         nh = 1
#     nw = int((wi/(4*h))*points[i])
#     print('points',points[i])
#     print('nt',nt)
#     print('nh',nh)
#     print('nw',nw)

    array, rel_dim = body.PA(xtrans, ytrans, farm)
    #array, rel_dim = body.attenuator(xtrans, ytrans, farm, D,nr,ntheta,nz)
    
    diff_result1, rad_result1, RAO_vals1, lam = solve.hydro(array, B, depth, w1, farm, controls)
    total1, incoming_fse1, x1, x2, nx, y1, y2, ny = solve.elevation(res, lam, diff_result1, rad_result1, RAO_vals1, farm, rad, controls, N, attenuator, rel_dim)
    rao1 = abs(total1[40,0])
    RAO_07.append(rao1)
    # rao1 = RAO_vals1[0, 0] # np.mean(total)
    # RAO_07.append(rao1)

    diff_result2, rad_result2, RAO_vals2, lam = solve.hydro(array, B, depth, w2, farm, controls)
    total2, incoming_fse2, x1, x2, nx, y1, y2, ny = solve.elevation(res, lam, diff_result2, rad_result2, RAO_vals2, farm, rad, controls, N, attenuator, rel_dim)
    rao2 = abs(total2[40,0])
    RAO_125.append(rao2)
    # rao2 = RAO_vals2[0, 0] # np.mean(total)
    # RAO_125.append(rao2)

    diff_result3, rad_result3, RAO_vals3, lam = solve.hydro(array, B, depth, w3, farm, controls)
    total3, incoming_fse3, x1, x2, nx, y1, y2, ny = solve.elevation(res, lam, diff_result3, rad_result3, RAO_vals3, farm, rad, controls, N, attenuator, rel_dim)
    rao3 = abs(total3[40,0])
    RAO_13.append(rao3)
    # rao3 = RAO_vals3[0, 0] # np.mean(total)
    # RAO_13.append(rao3)
    
    # Append the radial mesh value
    panels.append(points)

percent1 = [0]  
for i in range(1, len(RAO_07)):
    diff1 = 100 * (RAO_07[i] - RAO_07[i - 1]) / RAO_07[i - 1]
    percent1.append(diff1)
print('percent1',percent1)
percent2 = [0]  
for i in range(1, len(RAO_125)):
    diff2 = 100 * (RAO_125[i] - RAO_125[i - 1]) / RAO_125[i - 1]
    percent2.append(diff2)
print('percent2',percent2)
percent3 = [0]  
for i in range(1, len(RAO_13)):
    diff3 = 100 * (RAO_13[i] - RAO_13[i - 1]) / RAO_13[i - 1]
    percent3.append(diff3)
print('percent3',percent3)

#Blue: #377eb8 Orange: #ff7f00 Green: #4daf4a Purple: #984ea3 Yellow: #ffff33 Cyan: #a65628
# Plot the percent difference
plt.figure(figsize=(12, 6))
plt.plot(points, percent1, marker='o',color='#377eb8',label='$\omega$=1.0 rad/s')
plt.plot(points, percent2, marker='o',color='#ff7f00',label='$\omega$=1.25 rad/s')
plt.plot(points, percent3, marker='o',color='#984ea3',label='$\omega$=1.3 rad/s')
#plt.scatter(5.07,0.0878,marker='*',color='black',s=400,zorder=10,label='Chosen AR')
plt.xlabel('Aspect Ratio',fontsize=20)
plt.ylabel('Percent Change Wave Elevation [%]',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize=20, markerscale=1)
plt.tight_layout()
print('yip')
plt.savefig('free_surface_conv.pdf')
print('eee')
#plt.show()