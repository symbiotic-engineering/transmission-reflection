import numpy as np
import matplotlib.pyplot as plt
import body
import solve

B = 0                                   # wave heading [rad]
depth = 1.37                            # water depth at basin [m]
amp = 0.03                              # amplitude from experiments [m]
T = np.array([1.55])    #1.0,1.20,1.25,1.33,1.39]) # wave period from exp [s] np.linspace(1.0,1.33,30)
w = (2*np.pi)/T                         # wave frequency [rad/s]
reactive = False                        # whether controls are engaged
point_abs = True                       # device being modeled (True = PA, False = OS)
res = 1                                 # grid resolution multiplier
g = 9.81

if point_abs:
    # RAO presented in m/m
    singleRAO_exp = np.array([0.873208629567549,0.944754952740845,1.00243356061469,0.970693867686127,0.879420254492409])
    B_difference = np.array([2.601353059523+11.0361509650696*1j,0.937927576643199+6.28643744791731*1j,0.1571712306061+1.39353317578805*1j,0.4569846084097+4.06849014436796*1j,1.4500030088802+14.7873526149951*1j])

else:
    # RAO converted and represented in rad/m
    k = w**2/g
    singleRAO_exp = (np.array([320.908542699394,376.425396443474,379.078710327915,379.053394198879,392.783132781226]))/1000#*np.pi/180)
    print('exp RAO',singleRAO_exp)
    B_difference = np.array([5.21423039863538+1.81357171210956*1j,2.46409818130205+0.684169279546558*1j,2.14502513133324+0.579566593860094*1j,1.76241880789527+0.457753852086662*1j,1.48200819928462+0.386856344516752*1j]) 

ex_amp = []
RAO_amp = []
friction = []
added_mass = []
damping = []
viscous = []
natural_period = []

meshedbody = body.initialize(point_abs)              # create meshed body

for i in range(np.size(w)):                              
    RAO, diff_result, rad_result, exc_force, added_m, damp, viscous_damp, nat_per = solve.hydro(meshedbody,B,depth,w[i],reactive,point_abs,B_difference[i])   # compute RAOs, excitation force
    ex_real = np.real(exc_force)                # take the real part of the complex excitation force
    ex_arg = np.angle(exc_force)                # compute the argument of the trigonometric functions
    ex_unit_amp = (ex_real/np.cos(ex_arg))      # divide real part by cos(arg) to obtain amplitude

    if point_abs == False:
        RAO = np.abs(RAO)*k[i]*180/np.pi/1000 # this gives you deg/mm          # rotation angle turned to "RAO", phi [deg] *(k*amp)

    ex_amp.append(ex_unit_amp*amp)
    RAO_amp.append(np.abs(RAO))
    added_mass.append(added_m)
    damping.append(damp)
    viscous.append(viscous_damp)
    natural_period.append(nat_per)
print('RAO amp',RAO_amp)
print('added mass',added_mass)
print('damping',damping)

# plot excitation force wrt wave frequency
#plt.figure(figsize=(8,6))
T_exp = np.array([1.0,1.20,1.25,1.33,1.39])

#plt.plot(T,added_mass/np.max(added_mass),label='Added Mass',color='#CC79A7',marker='*',markersize=12,linewidth=2,linestyle=(0, (3, 1, 1, 1, 1, 1)))
#plt.plot(T,damping/np.max(damping),label='Radiation Damping',color='#56B4E9',marker='s',markersize=8,linewidth=2)
#plt.plot(T,ex_amp/np.max(ex_amp),label='Excitation Force',color='#F0E442',marker='o',markersize=8,linewidth=2)
#plt.plot(T,viscous/np.max(viscous),label='Viscous Damping',color='#009E73',marker='>',markersize=8,linewidth=2)

#plt.plot(T_exp,singleRAO_exp,label='Experiments',color='#56B4E9',marker='o',markersize=12,linewidth=2)
#plt.plot(T,RAO_amp,label='Model',color='#E69F00',marker='^',markersize=12,linewidth=2,linestyle=(0, (3, 1, 1, 1, 1, 1)))

plt.plot(T_exp,np.abs(B_difference),label='Magnitude',color='#CC79A7',marker='*',markersize=14,linewidth=3)
plt.plot(T_exp,np.real(B_difference),label='Real Part',color='#56B4E9',marker='o',markersize=12,linewidth=2)
plt.plot(T_exp,np.imag(B_difference),label='Imaginary Part',color='#009E73',marker='s',markersize=8,linewidth=2,linestyle=(0, (3, 1, 1, 1, 1, 1)))

plt.xticks(ticks=T_exp,fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Period [s]',fontsize=20)
plt.ylabel('Emprical Damping [$kg-m/s^2$]',fontsize=18)
plt.grid()
plt.legend(fontsize=15, markerscale=1)  #,loc='lower left')
plt.tight_layout()
plt.savefig('test.pdf')
plt.clf()