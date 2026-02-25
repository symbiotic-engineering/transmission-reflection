import numpy as np
import matplotlib.pyplot as plt

# Define constants
z = 0               # evaluating at z=0 [m]
w = 0.571            # dominant frequency [rad/s]
H_base = 2.19        # baseline wave height [m]
rho = 1025          # density of seawater [kg/m^3]
g = 9.81            # gravitational constant [m/s^2]
k = w**2 / g        # wave number
h = -40             # water depth
L = -h              # length
C_d = 1.0           # drag coefficient
C_m = 2.0           # inertia coefficient
rho_steel = 7980    # density of 316 stainless steel [kg/m^3]

# parameters that could change
# XL turbine dimensions (range for monopiles at depths > 30 m)
d_XL = np.linspace(7,10,5)                # diameter of XL monopiles [m]
t_XL = np.linspace(0.070,0.110,5)         # corresponding thickness of XL monopiles [m]

R = d_XL / 2             # outer radius [m]
thickness = t_XL         # monopile wall thickness [m]
r = R - thickness        # inner radius [m]
d = d_XL
x = -d / 2               # edge of turbine setting x=0 to turbine center [m]

# Calculate properties
volume = np.pi * (R**2 - r**2) * L                       # volume of submerged monopile [m^3]
mass = rho_steel * volume                                # mass of submerged monopile [kg]
I = ((mass * (R**2 + r**2)) / 4) + ((mass * L**2) / 12)  # moment of inertia about the y-axis (pitch) [kg-m^2]
y = R                                                    # distance to centroidal axis

T = 2 * np.pi / w                                        # wave period [s]
H_values = np.linspace(0.900,0.999,3) * H_base           # reduced wave heights, percentages chosen as examples 
percent_reduction_H = (1 - H_values / H_base) * 100      # percent reduction in H

def calculate_damage(H,I,y,d):
    def calculate_bending_stress(t,I,y,d):
        u = ((H * w) / 2) * (np.cosh(k * (z + h)) / np.sinh(k * h)) * np.cos(k * x - w * t)         # flow velocity
        u_dot = -((H * w**2) / 2) * (np.cosh(k * (z + h)) / np.sinh(k * h)) * np.sin(k * x - w * t) # flow acceleration
        f = C_m * rho * (np.pi / 4) * d**2 * u_dot + C_d * d * 0.5 * rho * u * np.abs(u)            # force of wave per unit length [N/m]
        bending_moment = f * (L + H / 2)                                                            # L + H/2 is the moment arm [N-m] or [kg-m^2/s^2]
        sigma = (y / I) * bending_moment                                                            # bending stress (Pa) [m/s^2]
        return sigma

    # Evaluate bending stress over the time range t = [0, T]
    t_values = np.linspace(0, T, 100)
    sigma_values = np.array([calculate_bending_stress(t,I,y,d) for t in t_values])

    # Find maximum and minimum sigma to obtain the range
    sigma_max = np.max(sigma_values)                    # this value is the amplitude of the sinusoidal stress curve
    sigma_min = np.min(sigma_values)                    # this is also the amplitude, but the negative value
    sigma_range = sigma_max - sigma_min                 # when you take this difference, you are effectively doubling the amplitude to get the sine wave height

    # a_bar value taken from experiments
    # for N > 10^6, m = 5
    m = 5
    a_bar = np.exp(13.617) 
    N_f = a_bar * (sigma_range**(-m))

    n_cycle = (25 * 8760 * 60 * 60) / T     # number of cycles in turbine lifetime (25 yrs)
    D = n_cycle / N_f
    return D

# Calculate baseline damage at H_base
percent_diff = [[],[],[],[],[]]
for i in range(np.size(R)):
    D_baseloop = calculate_damage(H_base,I[i],y[i],d[i])
    percent_diffinter = []
    for H in H_values:
        D_newloop = calculate_damage(H,I[i],y[i],d[i])
        percent_diffloop = 100 * ((D_newloop - D_baseloop)/ D_baseloop)

        percent_diffinter.append(percent_diffloop)
    percent_diff[i] = percent_diffinter

print('percent diff',percent_diff)

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(percent_reduction_H, percent_diff[0], marker='o',color='#D55E00',markersize=10)
plt.plot(percent_reduction_H, percent_diff[1], marker='>',markersize=10)
plt.plot(percent_reduction_H, percent_diff[2], marker='+',markersize=10)
plt.plot(percent_reduction_H, percent_diff[3], marker='x',markersize=10)
plt.plot(percent_reduction_H, percent_diff[4], marker='*',markersize=10)
plt.xlabel('Reduction of Wave Height [%]',fontsize=20)
plt.ylabel('Reduction in Fatigue Damage [%]',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.savefig('damage_reduction.pdf')
