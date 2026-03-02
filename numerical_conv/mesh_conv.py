'''this is probably the least straightforward and least informative
script in this whole repository. The geometric ratios I used for
each body are writted in my notes instead of leaving them commented
out in this file. I'm not sure the best way to consolidate this.'''

import body
import solve
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import wave_height

xtrans = [0,0]
ytrans = [-50,50]
controls = False
point_absorber = True
B = 0
depth = 40
x_center = 0
g = 9.81
w = np.flip((2*np.pi)/np.array([5,6,7,8,9,10,11,12]))

# Initialize lists to store results
panels = []
RAOS = [[],[],[],[],[],[],[],[]]
total_height = [[],[],[],[],[],[],[],[]]
perc_dif = [[],[],[],[],[],[],[],[]]
points = np.linspace(0.05,5,12)

for j in range(len(w)):
    RAO_loop = []
    total_loop = []
    percent_dif = []
    for i in range(len(points)):
        # for free surface
        res = points[i]
        print('res',res)

        # # for OSWEC
        wi, th, h = 18, 1.905, 10.8
        nt = int((th/1.75)*points[i])
        if nt < 1:
            nt = 1
        nh = int((h/4.5)*points[i])
        if nh < 1:
            nh = 1
        nw = int((wi/4.5)*points[i])
        # print('points',points[i])
        # print('nw',nw)
        # print('nt',nt)
        # print('nh',nh)

        # for PA
        # r,l = 14.6/2, 4.445
        # nr = int((r/(1.0*l))*points[i])
        # ntheta = int((np.pi/1.0)*points[i])
        # nz = int((2.0*l/r)*points[i])
        # print('AR',points[i])
        # print('nr',nr)
        # print('ntheta',ntheta)
        # print('nz',nz)
        nr, ntheta, nz = 16, 32, 12 

        # array, rel_dim, char_dim, budal_limit = body.PA(xtrans, ytrans, farm,w1,nr, ntheta, nz)
        array, rel_dim, char_dim = body.initialize(xtrans,ytrans,w[j],x_center,point_absorber,nr,ntheta,nz,nt,nw,nh)
    
        diff_result,rad_result,RAO_vals,lam,CWR = solve.hydro(array,B,depth,w[j],char_dim,controls,point_absorber)
        total, incoming_fse, x1, x2, nx, y1, y2, ny, elevation_of_interest = solve.elevation(res, lam, diff_result, rad_result, RAO_vals, controls, rel_dim)

        #ref1,trans1,EB1,KD1,power_abs1 = wave_height.wave_height(total1, incoming_fse1,xtrans,ytrans,farm,rel_dim,w1,nx,ny,x1,x2,y1,y2)
        k = w**2/g
        RAO_loop.append(RAO_vals[0])
        total_loop.append(elevation_of_interest)
        if i > 0:
            percent_dif_loop = (total_loop[i-1] - elevation_of_interest) / total_loop[i-1]
        else:
            percent_dif_loop = np.array([0])
        print(np.abs(percent_dif_loop))
        percent_dif.append(np.abs(percent_dif_loop))
    total_height[j] = total_loop
    RAOS[j] = RAO_loop
    perc_dif[j] = percent_dif
    print('perc dif',percent_dif)
    
    # Append the radial mesh value
    panels.append(points)

# Plot the percent difference
markerstyle = ['o','x','+','s','*','>','<','^']
colorblind = ['#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7', '#000000','#F5C200']
plt.figure(figsize=(12, 7))
for i in range(len(w)):
    plt.plot(points, np.abs(perc_dif[i]), marker=markerstyle[i],color=colorblind[i],label='$\omega$ = %.3f' %w[i])
plt.xlabel('Panel Number Increase Factor',fontsize=26)
plt.ylabel('Difference in Prediction',fontsize=26)
plt.legend(fontsize=16)#,loc='center right')
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.tight_layout()
print('yip')
plt.savefig('free_surface_res.pdf')
print('eee')
#plt.show()