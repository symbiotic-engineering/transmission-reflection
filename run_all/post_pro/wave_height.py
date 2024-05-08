import numpy as np
import matplotlib.pyplot as plt

xgrid = 3000                                # size of grid in x-direction
ygrid = 5000                                # size of grid in y-direction
mxc = 300                                   # number of grid points in x
myc = 300                                   # number of grid points in y
x_conversion = xgrid/(mxc + 1)
y_conversion = ygrid/(myc + 1)

Hs = np.loadtxt('/mnt/c/Users/ov162/transmission-reflection/data/wave_elevation.csv',delimiter=',')
Hs = Hs.reshape(mxc + 1,myc + 1)

x = [1490, 1590, 1690, 1420, 1520, 1620]    # x-position of bodies
ya = 4500                                   # y-position of row 1
yb = 4400                                   # y-position of row 2
H = 1.3832                                  # avg significant wave height [m]

x_investigated = int(1550/x_conversion)
y_investigated = int(4385/y_conversion)

midline_height = Hs[:y_investigated,x_investigated]
distance_from_array = np.linspace(0,y_investigated*y_conversion,np.size(midline_height))

waveheight_ezio = np.savetxt('waveheight_ezio.csv',midline_height,delimiter=',')

percent_decrease = ((H-midline_height)/H)*100
plt.plot(distance_from_array,percent_decrease)
plt.show()

plt.plot(distance_from_array,midline_height)
plt.show()

