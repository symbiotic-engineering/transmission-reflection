'''
This script uses data from The National Lab of the Rockies, NLR 
(previously the National Renewable Energy Lab, NREL) wind energy report from 2024
The LEVELIZED COST OF ENERGY (LCOE) of an individual offshore fixed bottom
wind turbine is estimated here. The script allows for changes in different 
parameters related to wave-induced fatigue damage reduction, namely:

1. for capital cost (CAPEX) parameters:
tower
substructure
turbine_installation
substruct_found_install
commissioning
decommissioning

2. for operations and maintenance cost (OPEX) parameters:
maintenance

'''
import numpy as np
import matplotlib.pyplot as plt

#### NREL EXAMPLE FARM CONTAINS 50, 12MW TURBINES. WE NORMALIZE BY THIS AND
# THEN CAN COMPUTE FOR THE SOUTHFORK CASE #####
nlr_numturbines = 50
nlr_rating = 12*1000                        # [kW]

####### CAPEX ##############
##### all parameters given in 2023 USD/kWh ######
###### turbine costs (for total farm) ######
rotor_nacelle = 1487
tower = 283
turbine_cost = (rotor_nacelle + tower)

##### balance of system (BOS) costs (per kWh!)#####
development = 121
project_management = 2
substructure = 232                          # this is the monopile
foundation = 556
array_cable_sys = 477
export_cable_sys = 532
grid_connection = 258
turbine_installation = 112
substruct_found_install = 172

BOS = development + project_management + substructure + foundation + array_cable_sys + export_cable_sys + grid_connection + turbine_installation + substruct_found_install

###### soft costs ########
construction_insurance = 55
decommissioning = 145
construction_financing = 240
procurement_contigency = 228
install_contingency = 289
commissioning = 55

soft_cost = construction_insurance + decommissioning + construction_financing + procurement_contigency + install_contingency + commissioning

########### fixed charge rate (FCR) #############
## this parameter annualizes the CAPEX, so proper units are actually: year^-1 
FCR = 0.0676                                    # nominal = 0.0874

CAPEX_ann = (turbine_cost + BOS + soft_cost) * FCR  # annualized, [$/kWh/year]
print('farm annual CAPEX, $/kWh/yr',CAPEX_ann)

########### OPEX ####################
######## all units given in $/kWh/year ##########
operations = 22                                 # pretax
maintenance = 113                               
design_life = 25                                # [years]

OPEX = operations + maintenance
print('farm annual OPEX, $/kWh/yr',OPEX)

##### annual energy production (AEP) #####
## this is given in units MWh/MW/yr (MWh-yr/MW)
farm_AEP = 4295                                 # [MWh-yr/MW]

#### levelized cost of energy (LCOE) #######
LCOE = ((CAPEX_ann + OPEX)*1000) / farm_AEP         # $/MWh
print('LCOE [$/MWh]',LCOE)

minLCOE = 100                                   # value if lifetime is 30 yrs
maxLCOE = 150                                   # value if lifetime is 15 yrs