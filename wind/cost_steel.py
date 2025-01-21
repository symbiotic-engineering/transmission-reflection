import numpy as np
import matplotlib.pyplot as plt

# cost breakdown for OWTs: https://www.sciencedirect.com/science/article/pii/S1364032118307342#bib61
# parameters
rho_steel = 7500        # density of steel [kg/m^3]

# South Fork turbine dimensions
d_SF = 10.97            # diameter of SF wind turbines [m]
t_SF = 0.125            # thickness of SF wind turbines [m]

# XL turbine dimensions (range for monopiles at depths > 30 m)
d_XLmin = np.linspace(7,10,5)             # diameter of smallest XL monopile [m]
t_XLmin = np.linspace(0.070,0.110,5)         # thickness of smallest XL monopile [m]

hub_height = 110        # hub height estimation from DOE [m] https://www.osti.gov/servlets/purl/1219141
depth = 40              # water depth at SF [m]
total_mp_length = hub_height + depth        # total monopile length [m]

# mass and volume of south fork monopiles
volume_SF = (((d_SF/2)**2 - ((d_SF - t_SF)/2)**2) * np.pi * total_mp_length)   # SF steel volume [m^3]
mass_SF = volume_SF * rho_steel       # SF monopile mass [kg]

volume_XL = []
mass_XL = []
# mass and volume of XL monopiles
for i in range(np.size(d_XLmin)):
    volume_min = (((d_XLmin[i]/2)**2 - ((d_XLmin[i] - t_XLmin[i])/2)**2) * np.pi * total_mp_length)
    min_mass = volume_min * rho_steel       # XL min monopile mass [kg]
    volume_XL.append(volume_min)
    mass_XL.append(min_mass)

# btwn 66% and 79% of the turbine's mass is steel https://escholarship.org/content/qt7p62d9rt/qt7p62d9rt.pdf
# cost of stainless steel fluctuates btwn $1500 and $6500 per ton throughout the year https://www.sciencedirect.com/science/article/pii/S0301420714000932
cost_of_steel =  (1500 + 6500)/2    # avg [$/ton]
tons_to_kg = 907.185        # conversion factor of US tons to kg
cost_in_kg = cost_of_steel/tons_to_kg   # cost of steel in [$/kg]

cost_of_SF = mass_SF * cost_in_kg
cost_of_XL = [mass * cost_in_kg for mass in mass_XL]

percent_diff = [((cost_of_SF - cost)/cost_of_SF)*100 for cost in cost_of_XL]
print('monopile diameter',d_XLmin)
print('monopile thickness', t_XLmin)
print('percent reduction in cost', percent_diff)

plt.figure(figsize=(8, 7))
plt.plot(d_XLmin,percent_diff)
plt.gca().invert_xaxis()
plt.ylabel('Reduction in Cost [%]',fontsize=20)
plt.xlabel('Monopile Diameter [m]',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('cost_reduction.pdf')