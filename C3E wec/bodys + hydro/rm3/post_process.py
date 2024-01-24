import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# analyzing wave elevation for rm3 
# housekeeping
k = 0.1117            #0.0529          # finite depth, k = 0.1117 in infinite depth
w = 1.047
#w = np.array([0.5,0.75,0.85,0.9,1.047,1.2,1.35,1.5,1.75,2.0])
#k = np.array([0.0254842,0.05733945,0.07364934,0.08256881,0.11174404,0.14678899,0.18577982,0.2293578,0.31218145,0.4077472])
T = (2*np.pi)/w       # period [s]
A_I = 1               # incoming wave amplitude
g = 9.81
z = 0
xR = np.linspace(-100,0,75)
xD = np.linspace(0,100,75)
xt = np.linspace(-100,100,150)
t = 0

# wave elevation from capytaine
capy_elev = "elevation_data_rm3.csv"
el = pd.read_csv(capy_elev, header=None)
zeta = el.values
z_up = zeta[74,:75]
z_down = zeta[74,75:]

# incident wave elevation
inc_el = "incoming_fse_data_bw.csv"
inc_el = pd.read_csv(inc_el, header=None)
inc = inc_el.applymap(lambda x: complex(x.replace(' ', '')))
zinc_up = -1*inc.iloc[:75,74]
zinc_down = -1*inc.iloc[75:,74]
Ainc_up = zinc_up*np.exp(1j*k*xR - 1j*w*t)
Ainc_down = zinc_down*np.exp(1j*k*xD - 1j*w*t)

# Finding the wave height --> reflection
xr_lt = -100
xr_ut = -25

# Find the indices where xR is between the specified thresholds
indices = np.where((xR >= xr_lt) & (xR <= xr_ut))[0]

# Check if there are any points within the specified range
if len(indices) > 0:
    # Find the index of the maximum and minimum values within the specified range
    max_index_zu = indices[np.argmax(z_up[indices])]
    min_index_zu = indices[np.argmin(z_up[indices])]

    # Calculate the maximum and minimum values
    max_zu = z_up[max_index_zu]
    min_zu = z_up[min_index_zu]

    # Calculate the wave height and reflection coefficient
    hu = max_zu - min_zu
    kr = (hu / (2 * A_I)) - 1

    print('Reflected wave height:', hu)
    print('Reflection Coefficient:', kr)
else:
    print("No points within the specified range.")

# transmission
xd_thresh = 25
zd_index = np.where(xD >= xd_thresh)[0][0]
max_index_zd = zd_index + np.argmax(z_down[zd_index:])
max_zd = z_down[max_index_zd]
min_index_zd = zd_index + np.argmin(z_down[zd_index:])
min_zd = z_down[min_index_zd]
hd = max_zd - min_zd
kt = hd/(2*A_I)
print('Transmitted wave height',hd)
print('Transmission Coefficient',kt)
y0 = np.zeros_like(xt)
plt.plot(xR,z_up,color='red',label='Upstream')
plt.plot(xD,z_down,color='black',label='Downstream')
plt.plot(xt,E_inc,color='blue',linestyle=':',label='Incident')
plt.plot(xt,y0,color='green',linestyle='--',label='Still Water Line')
plt.plot(xR[max_index_zu], max_zu, 'go',markersize=6) 
plt.plot(xR[min_index_zu], min_zu, 'go',markersize=6)  
plt.plot(xD[max_index_zd], max_zd, 'go',markersize=6)  
plt.plot(xD[min_index_zd], min_zd, 'go',markersize=6) 
plt.legend()
plt.xlabel('x [m]')
plt.ylabel('Wave Elevation [m]')
plt.savefig('rm3_elev.pdf')
plt.show()

#print('zd',z_down)
#print('xd',xD)
