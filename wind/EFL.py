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

# Convert each EFL value to MPa and store the results
EFL_ALL_MPa = [EFL for EFL in EFL_values]

# Print the results
print("EFL in unit MPa = ", EFL_ALL_MPa)

## Convert from 5s to 20 years

conversion = 12 * 60 * 24 * 30 *12 *20
conversion

twenty_EFL_ALL_MPa = [EFL_ALL_MPa * conversion for EFL_ALL_MPa in EFL_ALL_MPa]
twenty_EFL_ALL_MPa

# MPa, 20 year
first_last_difference = round(twenty_EFL_ALL_MPa[0] - twenty_EFL_ALL_MPa[-1], 2)
print("The difference of MPa for highest and lowest wave height is", first_last_difference)

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