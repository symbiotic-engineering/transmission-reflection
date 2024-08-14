import numpy as np
from scipy.integrate import dblquad

# Define the integrand for numerical integration
def integrand(t, z, w, x):
    H = 0.8
    rho = 1025
    g = 9.81  # gravitational constant
    k = w**2 / g
    h = -50
    u = ((H*w)/2) * (np.cosh(k*(z + h)) / np.sinh(k*h)) * np.cos(k*x - w*t)
    u_dot = -((H*w**2)/2) * (np.cosh(k*(z + h)) / np.sinh(k*h)) * np.sin(k*x - w*t)
    C_d = 1.0
    C_m = 2.0
    d = 10.97
    f = C_m*rho*(np.pi/4)*d**2 * u_dot + C_d*d*0.5*rho*u*np.abs(u)
    return f

# Define the limits for z and t
def z_limits(t):
    return [-50, 0.4]

def t_limits(w):
    return [0, 2*np.pi/w]  # T value for each w

# Sweep through a range of w values
w_values = np.array([0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.3])  # Example range of w values
x_value = 10.97
results = []

for w_value in w_values:
    T_value = (2*np.pi) / w_value
    result, error = dblquad(lambda t, z: integrand(t, z, w_value, x_value), 
                           *t_limits(w_value), lambda t: z_limits(t)[0], lambda t: z_limits(t)[1],
                           epsabs=1e-5, epsrel=1e-5)
    result_divided = result / T_value
    results.append((w_value, result_divided))

# Print or plot the results
for w_value, result_divided in results:
    print(f"w = {w_value:.2f}, Time-averaged force / T = {result_divided:.4f}")

# Optional: Plotting the results
import matplotlib.pyplot as plt

w_vals, forces = zip(*results)
plt.plot(w_vals, forces, marker='o')
plt.xlabel('Wave Frequency (w)')
plt.ylabel('Time-Averaged Force / T')
plt.title('Time-Averaged Force Divided by T vs. Wave Frequency')
plt.grid(True)
plt.show()