import numpy as np
import matplotlib.pyplot as plt

# for WEC costs
inflation_2014_to_2022 = 1.3

# Standalone WEC (From Neary et al. 2014)
WEC_units = np.array([1, 10, 50, 100])
WEC_ICC = np.array([61140, 21400, 14490, 13630]) *inflation_2014_to_2022     # [2023 $/kW] CAPEX
WEC_AOE = np.array([4070, 1150, 460, 330]) *inflation_2014_to_2022           # [2023 $/kW-yr] OPEX

plt.xlabel("WEC Units") 
plt.ylabel("CAPEX [million $/MW]") 
plt.plot(WEC_units, WEC_ICC, color ="red") 
plt.show()

plt.xlabel("WEC Units") 
plt.ylabel("OPEX [$/kW-yr]") 
plt.plot(WEC_units, WEC_AOE, color ="blue") 
plt.show()