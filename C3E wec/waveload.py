import numpy as np
import matplotlib.pyplot as plt
import math

# # parameters
# c_M = 3                 # inertia coefficient
# rho = 1025              # water density[kg/m^3]
# g = 9.81                # gravitational constant [m/s^2]
# D = 1                   # cylinder diameter
# H = 5                   # wave height
# w = 1.047               # wave frequency [rad/s]
# d = 40                  # water depth [m]
# c_D = 1                 # resistance coefficient
# t = np.linspace(0,300)  # time

# # dispersion relation
# k = g/w**2              # simplifying to deep water for now
# print(k)

# # coefficients
# C_F1 = ((-1/8)*c_M*rho*g*np.pi*(D**2)*H*np.sinh(k*d))/np.cosh(k*d)
# C_FD = ((1/16)*c_D*rho*g*D*(H**2)*2*k*d + np.sinh(2*k*d))/np.sinh(2*k*d)

# # force
# F = C_F1*np.sin(w*t) + C_FD*abs(np.cos(w*t))

H = np.linspace(0,10,num=30)

#H = np.array([0.9,1]) # 25% reduction in wave height
F = H + H**2
dF = 1 + 2*H
#print(F)            # 40% reduction in force


plt.xlabel('Wave Height [m]',fontsize=17)
plt.ylabel('Force [-]',fontsize=17)
plt.plot(H,F,color = 'red')
#plt.plot(H,dF,color='blue')
plt.savefig('H_vs_F.pdf')
plt.show()