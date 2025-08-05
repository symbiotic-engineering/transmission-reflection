import numpy as np
import matplotlib.pyplot as plt
import bodies
import solve

B = 0                                   # wave heading [rad]
depth = 1.37                            # water depth at basin [m]
input_amplitude = np.array([0.01,0.02,0.03,0.05])                  # amplitude from experiments [m]
T = np.array([1.0,1.20,1.25,1.33,1.39]) # wave period from exp [s]
w = (2*np.pi)/T                         # wave frequency [rad/s]
reactive = False                        # whether controls are engaged
res = 1                                 # grid resolution multiplier

array = bodies.initialize()             # create meshed array

ex_forceOS1 = []
ex_forceOS2 = []
ex_forcePA3 = []
ex_forcePA4 = []
   
for omega in w:                                       
    RAO, diff_result, rad_result, exc_force = solve.hydro(array,B,depth,omega,reactive)   # compute RAOs, excitation force
    ex_real = np.real(exc_force)                # take the real part of the complex excitation force
    ex_arg = np.angle(exc_force)                # compute the argument of the trigonometric functions
    ex_unit_amp = (ex_real/np.cos(ex_arg))      # dive real part by cos(arg) to obtain amplitude
    
    ex_forceOS1.append(ex_unit_amp[0])
    ex_forceOS2.append(ex_unit_amp[1])
    ex_forcePA3.append(ex_unit_amp[2])
    ex_forcePA4.append(ex_unit_amp[3])

ex_ampOS1 = []
ex_ampOS2 = []
ex_ampPA3 = []
ex_ampPA4 = []
for amp in input_amplitude: 
    ex_OS1 = []
    ex_OS2 = []
    ex_PA3 = []
    ex_PA4 = []
    for i in range(np.size(ex_forceOS1)):
        ex_OS1.append(ex_forceOS1[i]*amp)
        ex_OS2.append(ex_forceOS2[i]*amp)
        ex_PA3.append(ex_forcePA3[i]*amp)
        ex_PA4.append(ex_forcePA4[i]*amp)
    ex_ampOS1.append(ex_OS1)
    ex_ampOS2.append(ex_OS2)
    ex_ampPA3.append(ex_PA3)
    ex_ampPA4.append(ex_PA4)
    
print('OS1 excitation forces amplitude',ex_ampOS1)
print('OS2 excitation forces amplitude',ex_ampOS2)
print('PA3 excitation forces amplitude',ex_ampPA3)
print('PA4 excitation forces amplitude',ex_ampPA4)


# plot excitation force wrt wave frequency
for i in range(np.size(input_amplitude)):
    plt.plot(T,ex_ampOS1[i])
plt.savefig('exforce_OS1.pdf')
plt.clf()

for i in range(np.size(input_amplitude)):
    plt.plot(T,ex_ampOS2[i])
plt.savefig('exforce_OS2.pdf')
plt.clf()

for i in range(np.size(input_amplitude)):
    plt.plot(T,ex_ampPA3[i])
plt.savefig('exforce_PA3.pdf')
plt.clf()

for i in range(np.size(input_amplitude)):
    plt.plot(T,ex_ampPA4[i])
plt.savefig('exforce_PA4.pdf')
plt.clf()

# # compute wave elevation
# total, incoming_fse, grid, radiation, diffraction = solve.elevation(res,diff_result,rad_result,RAO)

# # plot wave elevation
# Z = np.abs(np.real(total)/np.real(incoming_fse))
# X = grid[0]
# Y = grid[1]
# pcm = plt.pcolormesh(X, Y, Z)
# plt.xlabel("x")
# plt.ylabel("y")
# colorbar = plt.colorbar()
# colorbar.set_label(r"Total Wave Elevation, $\eta$")
# #pcm.set_clim([-2, 2])
# plt.tight_layout()
# print('tip')
# plt.savefig('abstotalfield.pdf')
# print('top')

# note TO OLIVIA: (np.abs(total)/np.abs(incoming_fse)) - np.abs(total/incoming_fse) = numerically zero
# note TO OLIVIA: when modeling scaled up vs scaled down verisons, the PAs produce the same RAO, 
# but the OSWEC RAO changes VERY SIGNIFICANTLY, from around 4 to aroun 0.08