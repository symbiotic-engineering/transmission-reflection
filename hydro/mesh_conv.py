'''this is probably the least straightforward and least informative
script in this whole repository. The geometric ratios I used for
each body are writted in my notes instead of leaving them commented
out in this file. I'm not sure the best way to consolidate this.'''

import body
import solve
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wave_height

xtrans = [0,0]
ytrans = [-50,50]
farm = False
controls = False
attenuator = True
B = 0
depth = 40
w1 = 0.7
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
points = np.linspace(2,10,17)
#array, rel_dim, char_dim, budal_limit = body.PA(xtrans, ytrans, farm,w1)
#diff_result1,rad_result1,RAO_vals1,lam,CWR = solve.hydro(array,B,depth,w1,char_dim,farm,controls,point_absorber=True)
#diff_result2, rad_result2, RAO_vals2, lam,CWR = solve.hydro(array, B, depth, w2,char_dim, farm, controls,point_absorber=True)
#diff_result3, rad_result3, RAO_vals3, lam,CWR = solve.hydro(array, B, depth, w3, char_dim,farm, controls,point_absorber=True)

for i in range(len(points)):
    # for free surface
    #res = points[i]
    #print('res',res)
    # # for OSWEC
    # wi, th, h = 18, 1, 16
    # nt = int((th/3)*points[i])
    # if nt < 1:
    #     nt = 1
    # nh = int((h/4)*points[i])
    # if nh < 1:
    #     nh = 1
    # nw = int((wi/4)*points[i])
    # print('points',points[i])
    # print('nt',nt)
    # print('nh',nh)
    # print('nw',nw)
    # for attenuator
    r, l = 2, 14
    nr = int((r/3)*points[i])
    ntheta = int((4*np.pi/5)*points[i])
    nz = int((l/6)*points[i])
    print('points',points[i])
    print('nr',nr)
    print('ntheta',ntheta)
    print('nz',nz)
    # for PA
    # r,l = 10.5, 6
    # nr = int((r/(1.0*l))*points[i])
    # ntheta = int((np.pi/1.0)*points[i])
    # nz = int((2.0*l/r)*points[i])
    # print('points',points[i])
    # print('nr',nr)
    # print('ntheta',ntheta)
    # print('nz',nz)


    # array, rel_dim, char_dim, budal_limit = body.PA(xtrans, ytrans, farm,w1,nr, ntheta, nz)
    array, rel_dim, char_dim, budal_limit = body.attenuator(xtrans, ytrans, farm, D,w1,nr,ntheta,nz)
    
    diff_result1,rad_result1,RAO_vals1,lam,CWR = solve.hydro(array,B,depth,w1,char_dim,farm,controls,point_absorber=False)
    #total1, incoming_fse1, x1, x2, nx, y1, y2, ny = solve.elevation(res, lam, diff_result1, rad_result1, RAO_vals1, farm, rad, controls, N, attenuator, rel_dim)
    #ref1,trans1,EB1,KD1,power_abs1 = wave_height.wave_height(total1, incoming_fse1,xtrans,ytrans,farm,rel_dim,w1,nx,ny,x1,x2,y1,y2)
    #rao1 = trans1[0]
    #print('RAO!',rao1)
    #RAO_07.append(rao1)
    rao1 = RAO_vals1[0, 0]
    RAO_07.append(rao1)

    diff_result2, rad_result2, RAO_vals2, lam,CWR = solve.hydro(array, B, depth, w2,char_dim, farm, controls,point_absorber=False)
    #total2, incoming_fse2, x1, x2, nx, y1, y2, ny = solve.elevation(res, lam, diff_result2, rad_result2, RAO_vals2, farm, rad, controls, N, attenuator, rel_dim)
    #ref2,trans2,EB2,KD2,power_abs2 = wave_height.wave_height(total2, incoming_fse2,xtrans,ytrans,farm,rel_dim,w2,nx,ny,x1,x2,y1,y2)
    #rao2 = trans2[0]
    #print('RAO!',rao2)
    #RAO_125.append(rao2)
    rao2 = RAO_vals2[0, 0] # np.mean(total)
    RAO_125.append(rao2)

    diff_result3, rad_result3, RAO_vals3, lam,CWR = solve.hydro(array, B, depth, w3, char_dim,farm, controls,point_absorber=False)
    #total3, incoming_fse3, x1, x2, nx, y1, y2, ny = solve.elevation(res, lam, diff_result3, rad_result3, RAO_vals3, farm, rad, controls, N, attenuator, rel_dim)
    #ref3,trans3,EB3,KD3,power_abs3 = wave_height.wave_height(total3, incoming_fse3,xtrans,ytrans,farm,rel_dim,w3,nx,ny,x1,x2,y1,y2)
    #rao3 = trans3[0]
    #print('RAO!',rao3)
    #RAO_13.append(rao3)
    rao3 = RAO_vals3[0, 0] # np.mean(total)
    RAO_13.append(rao3)
    
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
#plt.figure(figsize=(15, 15))
plt.figure(figsize=(12, 7))
plt.plot(points, percent1, marker='o',color='#377eb8',label='$\omega$=0.7 rad/s')
plt.plot(points, percent2, marker='o',color='#ff7f00',label='$\omega$=1.25 rad/s')
plt.plot(points, percent3, marker='o',color='#984ea3',label='$\omega$=1.3 rad/s')
plt.scatter(7,0.0,marker='*',color='black',s=400,zorder=10,label='Chosen AR')
plt.xlabel('Aspect Ratio',fontsize=30)
plt.ylabel('Percent Change in RAO [%]',fontsize=30)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.legend(fontsize=25, markerscale=1)
plt.tight_layout()
print('yip')
plt.savefig('free_surface_conv.pdf')
print('eee')
#plt.show()