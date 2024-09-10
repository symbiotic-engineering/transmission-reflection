import numpy as np
import matplotlib.pyplot as plt

# Data
gauge_x = np.array([6.79, 6.79, 6.79, 6.985, 6.595, 6.205, 5.815, 5.425, 6.79, 6.4, 
                    6.01, 5.62, 6.4, 10.335])
gauge_y = -1 * np.array([2.845, 3.095, 3.295, 4.2, 4.2, 4.2, 4.2, 4.2, 4.608, 4.608, 4.608, 
                         4.608, 5.87, 4.2])

WECx = np.array([6.79, 6.4, 6.01])
WECy = -1*np.array([4.2, 4.2, 4.2])

# Plotting
plt.figure(figsize=(12, 8))
plt.scatter(gauge_x, gauge_y, c='#009E73', marker='o', edgecolor='black', label='Gauges',s=400)
plt.scatter(WECx, WECy, c='#D55E00', marker='*', edgecolor='black', label='WECs', s=750)
plt.xlabel('X',fontsize=35)
plt.ylabel('Y',fontsize=35)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.legend(fontsize=30)
plt.tight_layout()
plt.grid(True)
plt.savefig('marker_locs.pdf')
