#!/usr/bin/env python
# coding: utf-8

# ## Morison Equation 
# From (GIDEON_ROAN_THESIS)

# In[3]:


import math
import numpy as np
import matplotlib.pyplot as plt

rho_water = 1023  # density of seawater (kg/m^3)
g = 9.81 # gravity (m/s^2)
C_d = 1.0  # Drag coefficient (dimensionless)
C_m = 2.0  # Inertia coefficient (dimensionless)
D = 14 # Cylinder diameter (m)
A_projected = np.pi * (D/2)**2  # Projected area (m^2)
V_submerged = A_projected * 6  # Submerged volume (m^3), assume 6 m submerge the water
H = 10  # Significant wave height (m)
wave_period = 5  # Wave period (s)
lambda_wave = 50 * H  # Wavelength (m)
z = 0.0  # Vertical coordinate (m)
h = 50  # Water depth (m)
x = 0.0  # Horizontal position at which force is calculated (m)
z = 0 # calcuated at water surface


# Morrison equation
def morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t):
    # Calculate wave frequency (omega)
    omega = (2 * np.pi) / wave_period
    
    # Calculate wave number (k)
    k = (2 * np.pi) / lambda_wave
    
    # Horizontal wave particle velocity (u) equation (3.18)
    u = ((H * omega * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.cos(k * x - omega * t)
    
    # Derivative of u with respect to time for the acceleration (du_dt) equation (3.20)
    du_dt = ((-H * omega**2 * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.sin(k * x - omega * t)
    
    # Morrison equation to calculate wave force (equation 3.17)
    F_wave = 0.5 * rho_water * C_d * np.abs(u) * u * A_projected + rho_water * C_m * V_submerged * du_dt
    
    return F_wave


# Create an array of time points
time_points = np.linspace(0, wave_period, 100)

# Calculate wave force for each time point
wave_forces = [morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t) for t in time_points]

# Plot the graph
plt.figure(figsize=(12, 6))
plt.plot(time_points, wave_forces, label='Wave Force (N)')
plt.title('Wave Force over Time at Wave Height = 10 m')
plt.xlabel('Time (s)')
plt.ylabel('Wave Force (N)')
plt.legend()
plt.grid(True)
plt.show()


# In[4]:


print(wave_forces)


# In[5]:


force_at_each_second = {t: force for t, force in zip(np.round(time_points, 0), wave_forces) if t in np.arange(0, int(time_points[-1])+1)}

force_at_each_second


# In[6]:


wave_forces = np.array([10023221.571594613, 11682216.261883667, 13254378.474108752, 14734829.044353865, 16119999.722717851, 17407609.037098955, 18596618.44405161, 19687169.796646886, 20680505.454584427, 21578872.634511806, 22385413.8442451, 23104045.459430162, 23739326.681716517, 24296321.260872144, 24780454.467186477, 25197367.863369163, 25552774.44599678, 25852316.705111764, 26101430.08724517, 26305214.243002135, 26468314.297161117, 26594814.199371237, 26688143.99995514, 26751002.651551917, 26785297.667367216, 26787056.407815754, 26645730.794591997, 26317232.962854765, 25804575.23128724, 25112444.863333695, 24247142.988647424, 23216502.68643613, 22029787.470372036, 20697571.709740713, 19231604.79221348, 17644661.075767722, 15950377.887029655, 14163083.997376326, 12297621.143781062, 10369161.256506186, 8393022.108877875, 6384484.114726915, 4358610.966542874, 2330076.7325458415, 313001.9149869457, -1679199.1830402398, -3633957.626623231, -5539673.62333908, -7385823.227221989, -9163046.577377468, -10863216.577994147, -12479487.244543977, -14006321.270145576, -15439496.70109109, -16776092.946008027, -18014456.67365188, -19154148.474589866, -20195871.46690402, -21141383.310623243, -21993393.355337027, -22755446.876195073, -23431798.5515802, -24027277.49800776, -24547146.301685475, -24996956.569719292, -25382403.5658803, -25709182.495538022, -25982848.961881723, -26208686.031640835, -26391580.22459742, -26535908.57932111, -26645438.750428464, -26723243.863523085, -26771633.59657396, -26792102.675060995, -26739911.744023908, -26504729.96144305, -26083669.775690127, -25480587.270332627, -24700985.604571324, -23751943.379940137, -22642022.750521228, -21381158.66659903, -19980530.925097313, -18452420.957075253, -16810055.50894929, -15067439.56639243, -13239181.025019601, -11340309.727576291, -9386093.56161417, -7391854.3434324525, -5372786.202942177, -3343779.130304178, -1319250.2496249713, 687014.7507445188, 2662007.7754471144, 4593630.019824138, 6470807.401518792, 8283588.187723676, 10023221.571594605])

bins = np.arange(0, wave_period + 1, 1)  # Bins from 0 to T seconds
bin_centers = (bins[:-1] + bins[1:]) / 2  

bin_indices = np.digitize(time_points, bins, right=True)  # Assign time points to bins

# Initialize an array to hold the average force for each bin
average_forces = np.zeros(len(bins)-1)

# Calculate average force for each bin
for i in range(1, len(bins)):
    indices = bin_indices == i  
    if np.any(indices):
        average_forces[i-1] = np.mean(wave_forces[indices])  # Average force in the bin

average_force_per_bin = dict(zip(bin_centers, average_forces))

average_force_per_bin


# In[7]:


average_force_per_interval = {f"{i}-{i+1}": force for i, force in enumerate(average_forces)}

average_force_per_interval


# In[8]:


T_values = np.arange(1, 6)  # Denominator cannot be zero
omega_values = 2 * np.pi / T_values

force_values = list(average_force_per_interval.values())

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(omega_values, force_values, marker='o', linestyle='-', color='b')
plt.xlabel('Angular Frequency $\omega$ (rad/s)')
plt.ylabel('Force per Unit Length (N/m)')
plt.title('Force per Unit Length vs Angular Frequency For Wave height 10 m')
plt.grid(True)
plt.show()


# In[ ]:





# In[ ]:





# In[41]:


import numpy as np
import matplotlib.pyplot as plt


rho_water = 1023  # Density of seawater (kg/m^3)
g = 9.81  # Gravity (m/s^2)
C_d = 1.0  # Drag coefficient (dimensionless)
C_m = 2.0  # Inertia coefficient (dimensionless)
D = 14  # Cylinder diameter (m)
A_projected = np.pi * (D/2)**2  # Projected area (m^2)
submerge_depth = 6  # Assume 6 m below the water surface
V_submerged = A_projected * submerge_depth  # Submerged volume (m^3)
H = np.array([10, 9.5, 9, 8.5, 8, 7.5, 7]) # Array of significant wave heights (m)
wave_period = 5  # Wave period (s)
h = 50  # Water depth (m)
x = 0.0  # Horizontal position at which force is calculated (m)
z = 0.0  # Vertical coordinate at the water surface

# Morrison equation
def morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t):
    omega = (2 * np.pi) / wave_period
    
    # Calculate wave number (k)
    k = (2 * np.pi) / lambda_wave
    
    # Horizontal wave particle velocity (u) equation (3.18)
    u = ((H * omega * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.cos(k * x - omega * t)
    
    # Derivative of u with respect to time for the acceleration (du_dt) equation (3.20)
    du_dt = ((-H * omega**2 * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.sin(k * x - omega * t)
    
    # Morrison equation to calculate wave force (equation 3.17)
    F_wave = 0.5 * rho_water * C_d * np.abs(u) * u * A_projected + rho_water * C_m * V_submerged * du_dt
    
    return F_wave

# Create an array of time points
time_points = np.linspace(0, wave_period, 100)

# Plot setup
plt.figure(figsize=(10, 6))

# Calculate and plot wave force for each wave height
for H in H:
    lambda_wave = 50 * H
    wave_forces = [morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t) for t in time_points]
    plt.plot(time_points, wave_forces, label=f'Height {H:.1f}')

# Graph styling
plt.title('Wave Force over Time for Different Wave Heights')
plt.xlabel('Time (s)')
plt.ylabel('Wave Force (N)')
plt.legend(title='Wave Heights (m)', loc='lower left')
plt.grid(True)

# Show the plot
plt.tight_layout()
plt.show()


# In[46]:


import numpy as np
import matplotlib.pyplot as plt

rho_water = 1023  # Density of seawater (kg/m^3)
g = 9.81  # Gravity (m/s^2)
C_d = 1.0  # Drag coefficient (dimensionless)
C_m = 2.0  # Inertia coefficient (dimensionless)
D = 14  # Cylinder diameter (m)
A_projected = np.pi * (D/2)**2  # Projected area (m^2)
submerge_depth = 6  # Assume 6 m below the water surface
V_submerged = A_projected * submerge_depth  # Submerged volume (m^3)
H = np.array([10, 9.5, 9, 8.5, 8, 7.5, 7]) # Array of significant wave heights (m)
wave_period = 5  # Wave period (s)
h = 50  # Water depth (m)
x = 0.0  # Horizontal position at which force is calculated (m)
z = 0.0  # Vertical coordinate at the water surface

# Morrison equation
def morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t):
    omega = (2 * np.pi) / wave_period
    k = (2 * np.pi) / lambda_wave
    u = ((H * omega * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.cos(k * x - omega * t)
    du_dt = ((-H * omega**2 * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.sin(k * x - omega * t)
    F_wave = 0.5 * rho_water * C_d * np.abs(u) * u * A_projected + rho_water * C_m * V_submerged * du_dt
    return F_wave

# Create an array of time points
time_points = np.linspace(0, wave_period, 100)

# Plot setup
plt.figure(figsize=(10, 6))

# Calculate and plot wave force for each wave height
average_forces = []
for wave_height in H:
    lambda_wave = 50 * H  
    wave_forces = [morison_equation(C_d, C_m, A_projected, V_submerged, wave_height, wave_period, lambda_wave, z, h, x, t) for t in time_points]
    average_force = np.mean(np.abs(wave_forces))
    average_forces.append(average_force)

# Correcting the plot command to use the correct H array
plt.plot(H, average_forces, marker='o', linestyle='-', color='b')
plt.xlabel('Wave Height (H) [m]')
plt.ylabel('Average Force per Unit Length (N/m)')
plt.title('Average Force per Unit Length vs Wave Height')
plt.grid(True)
plt.show()


# In[ ]:





# ## Olivia data Morison equation

# In[19]:


import numpy as np
import matplotlib.pyplot as plt


rho_water = 1023  # Density of seawater (kg/m^3)
g = 9.81  # Gravity (m/s^2)
C_d = 1.0  # Drag coefficient (dimensionless)
C_m = 2.0  # Inertia coefficient (dimensionless)
D = 14  # Cylinder diameter (m)
A_projected = np.pi * (D/2)**2  # Projected area (m^2)
submerge_depth = 6  # Assume 6 m below the water surface
V_submerged = A_projected * submerge_depth  # Submerged volume (m^3)
Hs = np.array([1.396, 1.395, 1.394, 1.393, 1.392, 1.391, 1.390, 1.389, 1.388, 1.387, 1.386, 1.385, 1.384, 1.383]) # Array of significant wave heights (m)
wave_period = 5  # Wave period (s)
h = 50  # Water depth (m)
x = 0.0  # Horizontal position at which force is calculated (m)
z = 0.0  # Vertical coordinate at the water surface

# Function to calculate wave force using the Morrison equation
def morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t):
    omega = (2 * np.pi) / wave_period
    
    # Calculate wave number (k)
    k = (2 * np.pi) / lambda_wave
    
    # Horizontal wave particle velocity (u) equation (3.18)
    u = ((H * omega * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.cos(k * x - omega * t)
    
    # Derivative of u with respect to time for the acceleration (du_dt) equation (3.20)
    du_dt = ((-H * omega**2 * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.sin(k * x - omega * t)
    
    # Morrison equation to calculate wave force (equation 3.17)
    F_wave = 0.5 * rho_water * C_d * np.abs(u) * u * A_projected + rho_water * C_m * V_submerged * du_dt
    
    return F_wave

# Create an array of time points
time_points = np.linspace(0, wave_period, 100)

# Plot setup
plt.figure(figsize=(10, 6))

# Calculate and plot wave force for each wave height
for H in Hs:
    lambda_wave = 50 * H  # Calculate wavelength as 50 * significant wave height
    wave_forces = [morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t) for t in time_points]
    plt.plot(time_points, wave_forces, label=f'Height {H:.3f}m')

# Graph
plt.title('Wave Force over Time for Different Wave Heights')
plt.xlabel('Time (s)')
plt.ylabel('Wave Force (N)')
plt.legend(title='Wave Heights (m)', loc='lower left')
plt.grid(True)
plt.tight_layout()
plt.show()


# In[22]:


import numpy as np
import matplotlib.pyplot as plt

rho_water = 1023  # Density of seawater (kg/m^3)
g = 9.81  # Gravity (m/s^2)
C_d = 1.0  # Drag coefficient (dimensionless)
C_m = 2.0  # Inertia coefficient (dimensionless)
D = 14  # Cylinder diameter (m)
A_projected = np.pi * (D/2)**2  # Projected area (m^2)
submerge_depth = 6  # Assume 6 m below the water surface
V_submerged = A_projected * submerge_depth  # Submerged volume (m^3)
Hs = np.array([1.396, 1.395, 1.394, 1.393, 1.392, 1.391, 1.390, 1.389, 1.388, 1.387, 1.386, 1.385, 1.384, 1.383]) # Array of significant wave heights (m)
wave_period = 5  # Wave period (s)
h = 50  # Water depth (m)
x = 0.0  # Horizontal position at which force is calculated (m)
z = 0.0  # Vertical coordinate at the water surface

# Morrison equation
def morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t):
    omega = (2 * np.pi) / wave_period
    
    # Calculate wave number (k)
    k = (2 * np.pi) / lambda_wave
    
    # Horizontal wave particle velocity (u) equation (3.18)
    u = ((H * omega * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.cos(k * x - omega * t)
    
    # Derivative of u with respect to time for the acceleration (du_dt) equation (3.20)
    du_dt = ((-H * omega**2 * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.sin(k * x - omega * t)
    
    # Morrison equation to calculate wave force (equation 3.17)
    F_wave = 0.5 * rho_water * C_d * np.abs(u) * u * A_projected + rho_water * C_m * V_submerged * du_dt
    
    return F_wave

# Create an array of time points
time_points = np.linspace(0, wave_period, 100)

# Calculate and plot wave force for each wave height
average_forces = []
for H in Hs:
    lambda_wave = 50 * H  # Calculate wavelength as 50 * significant wave height
    wave_forces = [morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t) for t in time_points]
    average_force = np.mean(np.abs(wave_forces))
    average_forces.append(average_force)

plt.figure(figsize=(10, 6))
plt.plot(Hs, average_forces, marker='o', linestyle='-', color='b')
plt.xlabel('Wave Height (H) [m]')
plt.ylabel('Average Force per Unit Length (N/m)')
plt.title('Average Force per Unit Length vs Wave Height')
plt.grid(True)
plt.show()


# In[ ]:





# ## EFL 

# In[28]:


import math
import numpy as np
import matplotlib.pyplot as plt

rho_water = 1023  # density of seawater (kg/m^3)
g = 9.81 # gravity (m/s^2)
C_d = 1.0  # Drag coefficient (dimensionless)
C_m = 2.0  # Inertia coefficient (dimensionless)
D = 14 # Cylinder diameter (m)
A_projected = np.pi * (D/2)**2  # Projected area (m^2)
V_submerged = A_projected * 6  # Submerged volume (m^3), assume 6 m submerge the water
H = 10  # Significant wave height (m)
wave_period = 5  # Wave period (s)
lambda_wave = 50 * H  # Wavelength (m)
z = 0.0  # Vertical coordinate (m)
h = 50  # Water depth (m)
x = 0.0  # Horizontal position at which force is calculated (m)
z = 0 # calcuated at water surface


# Function to calculate wave force using the Morrison equation
def morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t):
    # Calculate wave frequency (omega)
    omega = (2 * np.pi) / wave_period
    
    # Calculate wave number (k)
    k = (2 * np.pi) / lambda_wave
    
    # Horizontal wave particle velocity (u) equation (3.18)
    u = ((H * omega * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.cos(k * x - omega * t)
    
    # Derivative of u with respect to time for the acceleration (du_dt) equation (3.20)
    du_dt = ((-H * omega**2 * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.sin(k * x - omega * t)
    
    # Morrison equation to calculate wave force (equation 3.17)
    F_wave = 0.5 * rho_water * C_d * np.abs(u) * u * A_projected + rho_water * C_m * V_submerged * du_dt
    
    return F_wave


# Create an array of time points
time_points = np.linspace(0, wave_period, 100)

# Calculate wave force for each time point
wave_forces = [morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t) for t in time_points]

# Plot the graph
plt.figure(figsize=(12, 6))
plt.plot(time_points, wave_forces, label='Wave Force (N)')
plt.title('Wave Force over Time at Wave Height = 10 m')
plt.xlabel('Time (s)')
plt.ylabel('Wave Force (N)')
plt.legend()
plt.grid(True)
plt.show()


# In[29]:


def morison_equation_to_stress(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t, D):
    F_wave = morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t)
    
    # Assuming the force is distributed across the projected area (A_projected),
    stress = F_wave / A_projected
    
    return stress

# Calculate stress for each time point
stresses = [morison_equation_to_stress(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t, D) for t in time_points]

# Plot the stress graph
plt.figure(figsize=(12, 6))
plt.plot(time_points, stresses, label='Stress (Pa)')
plt.title('Stress over Time at Wave Height = 10 m')
plt.xlabel('Time (s)')
plt.ylabel('Stress (Pa)')
plt.legend()
plt.grid(True)
plt.show()


# In[30]:


max_stress = np.max(stresses)
max_stress


# In[31]:


m = 3
N = 10**6

EFL = ((max_stress**m)/N)**(1/m)
print("EFL in unit Pa = ", EFL)


# In[32]:


EFL_MPa = EFL/10**6
print("EFL in unit MPa = ", EFL_MPa)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[27]:


import numpy as np
import matplotlib.pyplot as plt

rho_water = 1023  # Density of seawater (kg/m^3)
g = 9.81  # Gravity (m/s^2)
C_d = 1.0  # Drag coefficient (dimensionless)
C_m = 2.0  # Inertia coefficient (dimensionless)
D = 14  # Cylinder diameter (m)
A_projected = np.pi * (D/2)**2  # Projected area (m^2)
submerge_depth = 6  # Assume 6 m below the water surface
V_submerged = A_projected * submerge_depth  # Submerged volume (m^3)
Hs = np.array([1.396, 1.395, 1.394, 1.393, 1.392, 1.391, 1.390, 1.389, 1.388, 1.387, 1.386, 1.385, 1.384, 1.383]) # Array of significant wave heights (m)
wave_period = 5  # Wave period (s)
h = 50  # Water depth (m)
x = 0.0  # Horizontal position at which force is calculated (m)
z = 0.0  # Vertical coordinate at the water surface

# Morrison equation
def morison_equation(C_d, C_m, A_projected, V_submerged, Hs, wave_period, lambda_wave, z, h, x, t):
    omega = (2 * np.pi) / wave_period
    
    # Calculate wave number (k)
    k = (2 * np.pi) / lambda_wave
    
    # Horizontal wave particle velocity (u) equation (3.18)
    u = ((Hs * omega * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.cos(k * x - omega * t)
    
    # Derivative of u with respect to time for the acceleration (du_dt) equation (3.20)
    du_dt = ((-Hs * omega**2 * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.sin(k * x - omega * t)
    
    # Morrison equation to calculate wave force (equation 3.17)
    F_wave = 0.5 * rho_water * C_d * np.abs(u) * u * A_projected + rho_water * C_m * V_submerged * du_dt
    
    return F_wave

# Create an array of time points
time_points = np.linspace(0, wave_period, 100)

# Plot setup
plt.figure(figsize=(10, 6))

# Calculate and plot wave force for each wave height
for H in Hs:
    lambda_wave = 50 * H 
    wave_forces = [morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t) for t in time_points]
    plt.plot(time_points, wave_forces, label=f'Height {H:.3f}')

# Graph 
plt.title('Wave Force over Time for Different Wave Heights')
plt.xlabel('Time (s)')
plt.ylabel('Wave Force (N)')
plt.legend(title='Wave Heights (m)', loc='lower left')
plt.grid(True)
plt.tight_layout()
plt.show()


# In[ ]:





# In[13]:


import numpy as np
import matplotlib.pyplot as plt

rho_water = 1023  # Density of seawater (kg/m^3)
g = 9.81  # Gravity (m/s^2)
C_d = 1.0  # Drag coefficient (dimensionless)
C_m = 2.0  # Inertia coefficient (dimensionless)
D = 14  # Cylinder diameter (m)
A_projected = np.pi * (D/2)**2  # Projected area (m^2)
submerge_depth = 6  # Assume 6 m below the water surface
V_submerged = A_projected * submerge_depth  # Submerged volume (m^3)
Hs = np.array([1.396, 1.395, 1.394, 1.393, 1.392, 1.391, 1.390, 1.389, 1.388, 1.387, 1.386, 1.385, 1.384, 1.383]) # Array of significant wave heights (m)
wave_period = 5  # Wave period (s)
h = 50  # Water depth (m)
x = 0.0  # Horizontal position at which force is calculated (m)
z = 0.0  # Vertical coordinate at the water surface

# Morrison equation
def morison_equation(C_d, C_m, A_projected, V_submerged, Hs, wave_period, lambda_wave, z, h, x, t):
    omega = (2 * np.pi) / wave_period
    
    # Calculate wave number (k)
    k = (2 * np.pi) / lambda_wave
    
    # Horizontal wave particle velocity (u) equation (3.18)
    u = ((Hs * omega * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.cos(k * x - omega * t)
    
    # Derivative of u with respect to time for the acceleration (du_dt) equation (3.20)
    du_dt = ((-Hs * omega**2 * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.sin(k * x - omega * t)
    
    # Morrison equation to calculate wave force (equation 3.17)
    F_wave = 0.5 * rho_water * C_d * np.abs(u) * u * A_projected + rho_water * C_m * V_submerged * du_dt
    
    return F_wave

# Function to calculate stress from wave force
def calculate_stress(F_wave, A_projected):
    return F_wave / A_projected

# Create an array of time points
time_points = np.linspace(0, wave_period, 100)

# Plot setup
plt.figure(figsize=(10, 6))

# Calculate and plot stress for each wave height
for H in Hs:
    lambda_wave = 50 * H 
    wave_forces = np.array([morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t) for t in time_points])
    stresses = calculate_stress(wave_forces, A_projected)
    plt.plot(time_points, stresses, label=f'Height {H:.3f} m')

# Graph
plt.title('Stress over Time for Different Wave Heights')
plt.xlabel('Time (s)')
plt.ylabel('Stress (Pa)')
plt.legend(title='Wave Heights', loc='upper right')
plt.grid(True)
plt.tight_layout()
plt.show()


# In[20]:


import numpy as np

m = 3
N = 10**6
max_stresses = []
EFL_values = []

# Calculate max stress for each wave height and EFL using the max stress
for H in Hs:
    wave_forces = np.array([morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t) for t in time_points])
    stresses = calculate_stress(wave_forces, A_projected)
    max_stress = np.max(stresses)
    max_stresses.append(max_stress)
    # Calculate the EFL for the max stress of this wave height
    EFL = (((max_stress ** m) / N) ** (1 / m))/10**6
    EFL_values.append(EFL)

# Print the results
for i, H in enumerate(Hs):
    print(f"Wave Height: {H}m, EFL: {EFL_values[i]:.4e} MPa")


# In[21]:


# Convert each EFL value to MPa and store the results
EFL_ALL_MPa = [EFL for EFL in EFL_values]

# Print the results
print("EFL in unit MPa = ", EFL_ALL_MPa)


# In[22]:


## Convert from 5s to 20 years

conversion = 12 * 60 * 24 * 30 *12 *20
conversion


# In[23]:


twenty_EFL_ALL_MPa = [EFL_ALL_MPa * conversion for EFL_ALL_MPa in EFL_ALL_MPa]
twenty_EFL_ALL_MPa


# In[24]:


# MPa, 20 year
first_last_difference = round(twenty_EFL_ALL_MPa[0] - twenty_EFL_ALL_MPa[-1], 2)
print("The difference of MPa for highest and lowest wave height is", first_last_difference)


# In[ ]:





# In[26]:


import matplotlib.pyplot as plt

# Wave heights as given
Hs = np.array([1.396, 1.395, 1.394, 1.393, 1.392, 1.391, 1.390, 1.389, 1.388, 1.387, 1.386, 1.385, 1.384, 1.383])

# Plotting the EFL values against wave heights
plt.figure(figsize=(10, 6))
plt.plot(Hs, twenty_EFL_ALL_MPa, marker='o', linestyle='-', color='b')
plt.title('Equivalent Fatigue Load (EFL) (20 years) vs. Wave Height')
plt.xlabel('Wave Height (m)')
plt.ylabel('EFL (MPa)')
plt.grid(True)
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:




