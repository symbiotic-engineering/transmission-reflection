import numpy as np
import matplotlib.pyplot as plt
import body
import solve

B = 0                                   # wave heading [rad]
depth = 1.37                            # water depth at basin [m]
amp = 0.03                              # amplitude from experiments [m]
T = np.linspace(1.0,1.39,30) #np.array([1.0,1.20,1.25,1.33,1.39]) # wave period from exp [s] np.linspace(1.0,1.33,30)
w = (2*np.pi)/T                         # wave frequency [rad/s]
reactive = False                        # whether controls are engaged
point_abs = False                               # device being modeled (True = PA, False = OS)
res = 1                                 # grid resolution multiplier
g = 9.81

body = body.initialize(point_abs)              # create meshed body

ex_amp = []
RAO_amp = []
friction = []

for omega in w:                                       
    RAO, diff_result, rad_result, exc_force = solve.hydro(body,B,depth,omega,reactive,point_abs)   # compute RAOs, excitation force
    ex_real = np.real(exc_force)                # take the real part of the complex excitation force
    ex_arg = np.angle(exc_force)                # compute the argument of the trigonometric functions
    ex_unit_amp = (ex_real/np.cos(ex_arg))      # divide real part by cos(arg) to obtain amplitude
    k = (omega**2)/g 

    if point_abs == False:
        RAO = np.abs(RAO)*(180/np.pi)*amp             # rotation angle turned to "RAO", phi [deg] *(k*amp)

    ex_amp.append(ex_unit_amp[0]*amp)
    RAO_amp.append(np.abs(RAO[0]))

#print('rao amp',RAO_amp)
if point_abs:
    singleRAO_exp = np.array([0.0257692298239545,0.0274396442301012,0.0281018705744813,0.0297429038743760,0.0269462010272066])/amp
else:
    singleRAO_exp = np.array([9.36959308878599,10.7539296010859,10.8207487874797,11.4851088124678,11.9011123201484])

#percent_error = ((singleRAO_exp - RAO_amp)/singleRAO_exp)*100
#print('error',percent_error)

# ## for OSWEC, believe i need to use the coupled RAO from the surge AND pitch DOFs
# plot excitation force wrt wave frequency
T_exp = np.array([1.0,1.20,1.25,1.33,1.39])
plt.plot(T,RAO_amp,label='Model',color='#E69F00',marker='^',markersize=12,linewidth=2,linestyle=(0, (3, 1, 1, 1, 1, 1)))
plt.plot(T_exp,singleRAO_exp,label='Experiments',color='#56B4E9',marker='o',markersize=12,linewidth=2)
#plt.xticks(ticks=T_exp,fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Period [s]',fontsize=20)
plt.ylabel('$\phi$',fontsize=20)
plt.grid()
plt.legend(fontsize=15, markerscale=1)
plt.tight_layout()
plt.savefig('OS_friction.pdf')
plt.clf()