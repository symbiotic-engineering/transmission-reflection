import numpy as np

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
volume_max = (((d_SF/2)**2 - ((d_SF - t_SF)/2)**2) * np.pi * total_mp_length)   # SF steel volume [m^3]
max_mass = volume_max * rho_steel       # SF monopile mass [kg]

# mass and volume of XL monopiles
for d_XLmin in d_XLmin:
    for t_XLmin in t_XLmin
        volume_min = (((d_XLmin/2)**2 - ((d_XLmin - t_XLmin)/2)**2) * np.pi * total_mp_length)
        min_mass = volume_min * rho_steel       # XL min monopile mass [kg]

# btwn 66% and 79% of the turbine's mass is steel https://escholarship.org/content/qt7p62d9rt/qt7p62d9rt.pdf
# cost of stainless steel fluctuates btwn $1500 and $6500 per ton throughout the year
cost_of_steel =  (1500 + 6500)/2    # avg [$/ton]
tons_to_kg = 907.185        # conversion factor of US tons to kg
cost_in_kg = cost_of_steel/tons_to_kg   # cost of steel in [$/kg]

cost_of_SF = max_mass * cost_in_kg
cost_of_min = min_mass * cost_in_kg

percent_diff = ((cost_of_SF - cost_of_min)/cost_of_SF)*100

print('cost of south fork', cost_of_SF)
print('cost of min XL',cost_of_min)
print('percent reduction in cost', percent_diff)
