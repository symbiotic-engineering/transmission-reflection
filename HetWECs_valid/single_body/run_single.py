import numpy as np
import matplotlib.pyplot as plt
import body
import solve

B = 0                                   # wave heading [rad]
depth = 1.37                            # water depth at basin [m]
amp = 0.03                              # amplitude from experiments [m]
T = np.array([1.0,1.20,1.25,1.33,1.39]) # wave period from exp [s] np.linspace(1.0,1.33,30)
w = (2*np.pi)/T                         # wave frequency [rad/s]
reactive = False                        # whether controls are engaged
point_abs = False                       # device being modeled (True = PA, False = OS)
res = 1                                 # grid resolution multiplier
g = 9.81

if point_abs:
    singleRAO_exp = np.array([0.0257692298239545,0.0274396442301012,0.0281018705744813,0.0297429038743760,0.0269462010272066])/amp
    B_difference = np.array([2.9026955256015+14.5324995808832*1j,1.3770262250851+11.2925280918833*1j,0.9779699085638+9.4650510415704*1j,0.2392082712787+4.32758046520889*1j,1.2257112432888+14.6075107154072*1j])

else:
    singleRAO_exp = np.array([9.36959308878599,10.7539296010859,10.8207487874797,11.4851088124678,11.9011123201484])
    B_difference = np.array([0.651822802027485+2.39469341142627*1j,0.511076965538471+1.58836712214393*1j,0.48949439086433+1.4565182465286*1j,0.434116387575848+1.28064520449171*1j,0.401400971258321+1.17037545676909*1j])
    #np.array([45.4596115013629+9.36593227053133*1j,31.2701499378457+2.56532862325432*1j,29.5026131719652+1.98966304374386*1j,25.7866171205714+1.3220086179136*1j,23.6565806904476+1.00464110184365*1j])

ex_amp = []
RAO_amp = []
friction = []
added_mass = []
damping = []
viscous = []

body = body.initialize(point_abs)              # create meshed body

for i in range(np.size(w)):                                      
    RAO, diff_result, rad_result, exc_force, added_m, damp, viscous_damp = solve.hydro(body,B,depth,w[i],reactive,point_abs,B_difference[i])   # compute RAOs, excitation force
    ex_real = np.real(exc_force)                # take the real part of the complex excitation force
    ex_arg = np.angle(exc_force)                # compute the argument of the trigonometric functions
    ex_unit_amp = (ex_real/np.cos(ex_arg))      # divide real part by cos(arg) to obtain amplitude

    if point_abs == False:
        RAO = np.abs(RAO)*(180/np.pi)*amp             # rotation angle turned to "RAO", phi [deg] *(k*amp)


    ex_amp.append(ex_unit_amp*amp)
    RAO_amp.append(np.abs(RAO))
    added_mass.append(added_m)
    damping.append(damp)
    viscous.append(viscous_damp)
print('RAO amp',RAO_amp)

# plot excitation force wrt wave frequency
plt.figure(figsize=(10,6))
T_exp = np.array([1.0,1.20,1.25,1.33,1.39])

#plt.plot(T,added_mass/np.max(added_mass),label='Added Mass',color='#CC79A7',marker='*',markersize=12,linewidth=2,linestyle=(0, (3, 1, 1, 1, 1, 1)))
#plt.plot(T,damping/np.max(damping),label='Radiation Damping',color='#56B4E9',marker='s',markersize=8,linewidth=2)
#plt.plot(T,ex_amp/np.max(ex_amp),label='Excitation Force',color='#F0E442',marker='o',markersize=8,linewidth=2)
#plt.plot(T,viscous/np.max(viscous),label='Viscous Damping',color='#009E73',marker='>',markersize=8,linewidth=2)

plt.plot(T_exp,singleRAO_exp,label='Experiments',color='#56B4E9',marker='o',markersize=12,linewidth=2)
plt.plot(T,RAO_amp,label='Model',color='#E69F00',marker='^',markersize=12,linewidth=2,linestyle=(0, (3, 1, 1, 1, 1, 1)))

#plt.plot(T_exp,np.abs(B_difference),label='Magnitude',color='#CC79A7',marker='*',markersize=14,linewidth=3)
#plt.plot(T_exp,np.real(B_difference),label='Real Part',color='#56B4E9',marker='o',markersize=12,linewidth=2)
#plt.plot(T_exp,np.imag(B_difference),label='Imaginary Part',color='#009E73',marker='s',markersize=8,linewidth=2,linestyle=(0, (3, 1, 1, 1, 1, 1)))

plt.xticks(ticks=T_exp,fontsize=15)
#plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Period [s]',fontsize=20)
plt.ylabel('Difference in Damping',fontsize=20)
plt.grid()
plt.legend(fontsize=15, markerscale=1)  #,loc='lower left')
plt.tight_layout()
plt.savefig('test.pdf')
plt.clf()