import numpy as np
import matplotlib.pyplot as plt

# Define constants
z = 0               # evaluating at z=0 [m]
d = 10.97           # diameter of turbine
x = -d / 2          # edge of turbine setting x=0 to turbine center [m]
w = 1.25            # dominant frequency [rad/s]
H_base = 0.8        # baseline wave height [m]
rho = 1025          # density of seawater [kg/m^3]
g = 9.81            # gravitational constant [m/s^2]
k = w**2 / g        # wave number
h = -40             # water depth
L = -h              # length
C_d = 1.0           # drag coefficient
C_m = 2.0           # inertia coefficient
R = d / 2           # outer radius [m]
thickness = 0.125   # monopile wall thickness [m]
r = R - thickness   # inner radius [m]
rho_steel = 7980    # density of 316 stainless steel [kg/m^3]

# Calculate properties
volume = np.pi * (R**2 - r**2) * L      # volume of submerged monopile [m^3]
mass = rho_steel * volume               # mass of submerged monopile [kg]
I = ((mass * (R**2 + r**2)) / 4) + ((mass * L**2) / 12)  # moment of inertia about the y-axis (pitch) [kg-m^2]
y = R                             # distance to centroidal axis

T = 2 * np.pi / w                 # period [s]
H_values = np.linspace(0.94, 0.99, 100) * H_base
percent_reduction_H = (1 - H_values / H_base) * 100  # percent reduction in H

def calculate_damage(H):
    def calculate_bending_stress(t):
        u = ((H * w) / 2) * (np.cosh(k * (z + h)) / np.sinh(k * h)) * np.cos(k * x - w * t)         # flow velocity
        u_dot = -((H * w**2) / 2) * (np.cosh(k * (z + h)) / np.sinh(k * h)) * np.sin(k * x - w * t) # flow acceleration
        f = C_m * rho * (np.pi / 4) * d**2 * u_dot + C_d * d * 0.5 * rho * u * np.abs(u)   # force of wave per unit length [N/m]
        bending_moment = f * (L + H / 2)                                   # L + H/2 is the moment arm [N-m] or [kg-m^2/s^2]
        sigma = (y / I) * bending_moment  # bending stress (Pa) [m/s^2]
        return sigma

    # Evaluate bending stress over the time range t = [0, T]
    t_values = np.linspace(0, T, 100)
    sigma_values = np.array([calculate_bending_stress(t) for t in t_values])

    # Find maximum and minimum sigma
    sigma_max = np.max(sigma_values)
    sigma_min = np.min(sigma_values)
    sigma_range = sigma_max - sigma_min

    # a_bar value taken from experiments
    # for N > 10^6, m = 5
    m = 5
    a_bar = np.exp(13.617) 
    N_f = a_bar * (sigma_range**(-m))

    n_cycle = (25 * 8760 * 60 * 60) / T     # number of cycles in turbine lifetime (25 yrs)
    D = n_cycle / N_f
    return D

# Calculate baseline damage at H_base
D_base = calculate_damage(H_base)

# Calculate percent difference of D for each H
percent_diff_D = [(calculate_damage(H) - D_base) / D_base * 100 for H in H_values]

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(percent_reduction_H, percent_diff_D, marker='o',color='#D55E00')
plt.xlabel('Reduction of Wave Height [%]',fontsize=20)
plt.ylabel('Reduction in Fatigue Damage [%]',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.savefig('damage_reduction.pdf')

# import numpy as np

# # Define constants
# z = 0               # evaluating at z=0 [m]
# d = 10.97           # diameter of turbine
# x = -d / 2          # edge of turbine setting x=0 to turbine center [m]
# w = 1.25            # dominant frequency [rad/s]
# H = 0.8             # wave height [m]
# rho = 1025          # density of seawater [kg/m^3]
# g = 9.81            # gravitational constant [m/s^2]
# k = w**2 / g        # wave number
# h = -40             # water depth
# L = -h              # length
# C_d = 1.0           # drag coeff
# C_m = 2.0           # inertia coeff
# R = d / 2           # outer radius [m]
# thickness = 0.125   # monopile wall thickness [m]
# r = R - thickness   # inner radius [m]
# rho_steel = 7980    # density of 316 stainless steel [kg/m^3]

# # Calculate properties
# volume = np.pi * (R**2 - r**2) * L      # volume of submerged monopile [m^3]
# mass = rho_steel * volume               # mass of submerged monopile [kg]
# I = ((mass * (R**2 + r**2)) / 4) + ((mass * L**2) / 12)  # moment of inertia about the y-axis (pitch) [kg-m^2]
# y = R                             # distance to centroidal axis

# T = 2 * np.pi / w                 # period [s]
# H_values = np.linspace(0.94,0.99,100)*0.8

# for H in H_values:
#     # Define the bending stress calculation function
#     def calculate_bending_stress(t):
#         u = ((H*w)/2) * (np.cosh(k*(z + h)) / np.sinh(k*h)) * np.cos(k*x - w*t)         # flow velocity
#         u_dot = -((H*w**2)/2) * (np.cosh(k*(z + h)) / np.sinh(k*h)) * np.sin(k*x - w*t) # flow acceleration
#         f = C_m*rho*(np.pi/4)*d**2 * u_dot + C_d*d*0.5*rho*u*np.abs(u)   # force of wave per unit length [N/m]
#         bending_moment = f * (L + H/2)                                   # L + H/2 is the moment arm [N-m] or [kg-m^2/s^2]
#         sigma = (y / I) * bending_moment  # bending stress (Pa) [m/s^2]
#         return sigma

#     # Evaluate bending stress over the time range t = [0, T]
#     t_values = np.linspace(0, T, 100)
#     sigma_values = np.array([calculate_bending_stress(t) for t in t_values])

#     # Find maximum and minimum sigma
#     sigma_max = np.max(sigma_values)
#     sigma_min = np.min(sigma_values)
#     sigma_range = sigma_max - sigma_min
#     #print(f"Bending Stress Range: {sigma_range:.4f} Pa")

#     # a_bar value taken from experiments
#     # for N > 10^6, m = 5
#     m = 5
#     a_bar = np.exp(13.617) 
#     N_f = a_bar*(sigma_range**(-m))

#     n_cycle = (25 * 8760 * 60 * 60) / T     # number of cycles in turbine lifetime (25 yrs)
#     D = n_cycle/N_f
#     print('total damage',D)