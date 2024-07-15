def limits(RAO_controlled,w,B_pto,char_dim,CWR,KD,power):
    import numpy as np
    
    # power available in wave
    rho = 1025                  # [kg/m^3] density of sea water
    g = 9.81                    # [m/s^2] gravitational constant
    amp = 1                     # [m] unit amplitude
    power_avail = (rho * g**2 * amp**2) / (4 * w)      # [kW/m]
    while np.any(CWR > abs(KD)):
        # Create a mask for elements that exceed the budal_limit
        mask = CWR > abs(KD)
        RAO_controlled = RAO_controlled.copy()
        # Update RAO_controlled for those specific elements
        RAO_controlled = 0.95 * RAO_controlled
        # Recalculate power for those specific elements
        new_power = 0.5 * B_pto * (abs(RAO_controlled * w * 1j))**2
        
        # Use the mask to update only the relevant elements
        power = np.where(mask, new_power, power)
        CWR = np.where(mask, power / (power_avail*char_dim), CWR)
    
    print('final capture width ratio',CWR)
    #print('final saturated power', power)

    return RAO_controlled, CWR