import numpy as np
import matplotlib.pyplot as plt
# Constants
xgrid = 3000
ygrid = 5000
mxc = 300
myc = 300
x_conversion = xgrid / (mxc + 1)
y_conversion = ygrid / (myc + 1)
x_positions = [1490, 1590, 1690, 1420, 1520, 1620]
ya = 4500
yb = 4400
H = 1.3832
x_investigated = int(1550 / x_conversion)
y_investigated = int(4385 / y_conversion)

# Constants
rho_water = 1023  # Density of seawater (kg/m^3)
g = 9.81  # Gravity (m/s^2)
C_d = 1.0  # Drag coefficient (dimensionless)
C_m = 2.0  # Inertia coefficient (dimensionless)
D = 10.97  # monopile diameter (m)
A_projected = np.pi * (D/2)**2  # Projected area (m^2)
submerge_depth = 40  # Assume 6 m below the water surface
V_submerged = A_projected * submerge_depth  # Submerged volume (m^3)
#Hs = np.loadtxt('/mnt/c/Users/ov162/transmission-reflection/data/OSWEC_elevation.csv', delimiter=',')
#Hs = Hs[:y_investigated, x_investigated]
#H_inc = np.loadtxt('/mnt/c/Users/ov162/transmission-reflection/data/blank_elevation.csv',delimiter=',')
#H_inc = H_inc[:y_investigated, x_investigated]

wave_period = 5  # Wave period (s)
h = 40  # Water depth (m)
x = 0.0  # Horizontal position at which force is calculated (m)
z = 0.0  # Vertical coordinate at the water surface
lambda_wave=0
# Morrison equation
def morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t):
    omega = (2 * np.pi) / wave_period
    
    # Calculate wave number (k)
    #k = (2 * np.pi) / lambda_wave
    k = omega**2/g
    
    # Horizontal wave particle velocity (u)
    u = ((H * omega * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.cos(k * x - omega * t)
    
    # Derivative of u with respect to time for the acceleration (du_dt)
    du_dt = ((-H * omega**2 * np.cosh(k * (z + h))) / (2 * np.sinh(k * h))) * np.sin(k * x - omega * t)
    
    # Morrison equation to calculate wave force
    F_wave = (0.5 * rho_water * C_d * np.abs(u) * u * A_projected) + (rho_water * C_m * V_submerged * du_dt)
    #other = (0.5 * rho_water * C_d * A_projected) + (rho_water * C_m * V_submerged)
    #print('other',other)
    #velocity = u*abs(u) + du_dt
    #F_wave=velocity
    return F_wave

# Function to calculate stress from wave force
def calculate_stress(F_wave, A_projected):
    return F_wave / A_projected

# Create an array of time points
time_points = np.linspace(0, wave_period, 100)
m = 3
N = 10**6
max_stresses = []
EFL_values = []

h_dummy = np.linspace(0,15,15)
wave_force_dummy =[]
for H in h_dummy:
    wave_forces_dummy = np.array([morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t) for t in time_points])
    #plt.plot(time_points,wave_forces_dummy)
    #plt.show()
    wave_force_dummy.append(max(wave_forces_dummy))
plt.plot(wave_force_dummy)
plt.show()

# Calculate max stress for each wave height and EFL using the max stress
for H in Hs:
    lambda_wave = 50 * H
    wave_forces = np.array([morison_equation(C_d, C_m, A_projected, V_submerged, H, wave_period, lambda_wave, z, h, x, t) for t in time_points])
    stresses = calculate_stress(wave_forces, A_projected)
    max_stress = np.max(stresses)
    max_stresses.append(max_stress)
    # Calculate the EFL for the max stress of this wave height
    EFL = (((max_stress ** m) / N) ** (1 / m))/10**6
    EFL_values.append(EFL)


max_stresses_i = []
EFL_values_i = []
# calculate stresses from incident wave
for H_i in H_inc:
    lambda_wave_i = 50 * H_i
    wave_forces_i = np.array([morison_equation(C_d, C_m, A_projected, V_submerged, H_i, wave_period, lambda_wave, z, h, x, t) for t in time_points])
    stresses_i = calculate_stress(wave_forces_i, A_projected)
    max_stress_i = np.max(stresses_i)
    max_stresses_i.append(max_stress_i)
    # Calculate the EFL for the max stress of this wave height
    EFL_i = (((max_stress_i ** m) / N) ** (1 / m))/10**6
    EFL_values_i.append(EFL_i)

# Convert from 5s to 20 years
conversion = 12 * 60 * 24 * 30 *12 *20
twenty_year_EFL = [EFL_values * conversion for EFL_values in EFL_values]
twenty_year_EFL_i =  [EFL_values_i * conversion for EFL_values_i in EFL_values_i]

# MPa, 20 year
percent_EFL = []
percent_Hs = []
for i in range(len(Hs)):
    first_last_difference = 100*((twenty_year_EFL_i[i] - twenty_year_EFL[i])/twenty_year_EFL_i[i])
    first_last_Hs_difference = 100*((H_inc[i] - Hs[i])/H_inc[i])
    percent_EFL.append(first_last_difference)
    percent_Hs.append(first_last_Hs_difference)
    #print("The percent difference of EFL for highest and lowest wave height is", first_last_difference)
    #print("The percent difference of the highest and lowest wave height is", first_last_Hs_difference)


# Plotting the EFL values against wave heights
distance_from_array = abs(np.linspace(0, y_investigated * y_conversion, len(Hs)) - ygrid)
plt.figure(figsize=(10, 6))
plt.plot(distance_from_array, percent_EFL, marker='o', linestyle='-', color='b')
plt.xlabel('Wave Height Reduction [%]')
plt.ylabel('EFL Reduction [%]')
plt.grid(True)
plt.savefig('/mnt/c/Users/ov162/transmission-reflection/wind/OS_perc_red_dist.pdf')
plt.show()