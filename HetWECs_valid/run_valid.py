import numpy as np
import matplotlib.pyplot as plt
import bodies
import solve

B = 0                                   # wave heading [rad]
depth = 1.37                            # water depth at basin [m]
T = np.array([1.0,1.25,1.33]) # wave period from exp [s]
w = (2*np.pi)/T                         # wave frequency [rad/s]
g = 9.81
k = w**2/g
reactive = False                        # whether controls are engaged
res = 1                                 # grid resolution multiplier

array = bodies.initialize()             # create meshed array

RAO_OS1 = []
RAO_OS2 = []
RAO_PA3 = []
RAO_PA4 = []

add_OS1 = []
add_OS2 = []
add_PA3 = []
add_PA4 = []

damp_OS1 = []
damp_OS2 = []
damp_PA3 = []
damp_PA4 = []

ex_OS1 = []
ex_OS2 = []
ex_PA3 = []
ex_PA4 = []

BPTO_OS1 = []
BPTO_OS2 = []
BPTO_PA3 = []
BPTO_PA4 = []

power_OS1 = []
power_OS2 = []
power_PA3 = []
power_PA4 = []

# empirically derived damping coeffs from isolated device tests
B_diffPA = np.array([2.601353059523+11.0361509650696*1j,0.1571712306061+1.39353317578805*1j,0.4569846084097+4.06849014436796*1j])
B_diffOS = np.array([5.21423039863538+1.81357171210956*1j,2.14502513133324+0.579566593860094*1j,1.76241880789527+0.457753852086662*1j])
# B_diffPA = np.array([2.601353059523+11.0361509650696*1j,0.937927576643199+6.28643744791731*1j,0.1571712306061+1.39353317578805*1j,0.4569846084097+4.06849014436796*1j,1.4500030088802+14.7873526149951*1j])
# B_diffOS = np.array([5.21423039863538+1.81357171210956*1j,2.46409818130205+0.684169279546558*1j,2.14502513133324+0.579566593860094*1j,1.76241880789527+0.457753852086662*1j,1.48200819928462+0.386856344516752*1j])

# empirically derived damping coeffs from array tests, PTO disengaged
Bd_vir = np.array([28.660523034874206+126.9728540816521*1j,27.20076324837215+81.03140354379195*1j,24.826928807253346+72.14867163750871*1j])
Bd_fran = np.array([23.291842420661137+135.57949865296885j,29.47379714707714+78.06298269229453*1j,25.50907937179154+69.8867808578375*1j])
Bd_meg = np.array([1.2266677151860503-93.80138359019993j,-81.22483790518561-124.79505322193616*1j,-104.14035743647037-121.02988347419092*1j])
Bd_jo = np.array([-18.564006005286327-145.8236468237688j,-117.2861077516863-127.67332609452676*1j,-132.97385798117375-116.30957769390076*1j])
# Bd_vir = np.array([28.660523034874206+126.9728540816521*1j,27.924842028513627+88.61001683278265*1j,27.20076324837215+81.03140354379195*1j,24.826928807253346+72.14867163750871*1j,17.463389920514537+67.42201135724427*1j])
# Bd_fran = np.array([23.291842420661137+135.57949865296885j,29.172103691844494+85.38048193975598*1j,29.47379714707714+78.06298269229453*1j,25.50907937179154+69.8867808578375*1j,16.83931066768374+65.94098934998279*1j])
# Bd_meg = np.array([1.2266677151860503-93.80138359019993j,-63.27161920377039-130.74541394657783*1j,-81.22483790518561-124.79505322193616*1j,-104.14035743647037-121.02988347419092*1j,-98.82567271173829-126.03085146322015*1j])
# Bd_jo = np.array([-18.564006005286327-145.8236468237688j,-82.16213144252684-126.85806877585861*1j,-117.2861077516863-127.67332609452676*1j,-132.97385798117375-116.30957769390076*1j,-108.10513921859607-117.18251379263745*1j])

# ## for PTO disengaged
# megRAOexp = np.array([0.819664657814324,0.916166493356739,0.887306516520389,0.863314218351526,0.992867758219807])
# joRAOexp = np.array([0.477702270828799,0.848824917334879,0.680824985210178,0.713990464053677,0.958269910240572])
# virRAOexp = np.array([323.381890238651,361.268452383110,346.366992147017,322.459657109623,385.006028894270])*(np.pi/180)/k
# franRAOexp = np.array([393.634896229757,384.644580468920,338.447782918374,316.329978665478,396.131175641520])*(np.pi/180)/k

## RAOS for PTO engaged
megRAOexp = np.array([0.773495530727839,0.770132236212850,0.700512165810569])
joRAOexp = np.array([0.258407324729989,0.106376997648463,0.0307451184177340])
virRAOexp = np.array([184.572145311627,195.391623861709,165.997190210489])*(np.pi/180)/k
franRAOexp = np.array([131.852630012738,153.262515441623,241.872472048376])*(np.pi/180)/k
   
for i in range(np.size(w)):    
    print('running wave period: ',T[i])                                   
    RAO, diff_result, rad_result, ex_force, added_mass, damping, B_PTOemp, mech_power = solve.hydro(array,B,depth,w[i],reactive,
                                                                            B_diffPA[i],B_diffOS[i],megRAOexp[i],
                                                                            joRAOexp[i],virRAOexp[i],franRAOexp[i],
                                                                            Bd_vir[i],Bd_fran[i],Bd_meg[i],Bd_jo[i])   # compute RAOs, excitation force
    
    RAO_OS1.append(np.abs(RAO[0])*k[i]*180/np.pi/1000)
    RAO_OS2.append(np.abs(RAO[1])*k[i]*180/np.pi/1000)
    RAO_PA3.append(np.abs(RAO[2]))
    RAO_PA4.append(np.abs(RAO[3]))

    add_OS1.append(added_mass[0])
    add_OS2.append(added_mass[1])
    add_PA3.append(added_mass[2])
    add_PA4.append(added_mass[3])

    damp_OS1.append(damping[0])
    damp_OS2.append(damping[1])
    damp_PA3.append(damping[2])
    damp_PA4.append(damping[3])

    ex_OS1.append(ex_force[0])
    ex_OS2.append(ex_force[1])
    ex_PA3.append(ex_force[2])
    ex_PA4.append(ex_force[3])

    BPTO_OS1.append(B_PTOemp[0])
    BPTO_OS2.append(B_PTOemp[1])
    BPTO_PA3.append(B_PTOemp[2])
    BPTO_PA4.append(B_PTOemp[3])

    power_OS1.append(mech_power[0])
    power_OS2.append(mech_power[1])
    power_PA3.append(mech_power[2])
    power_PA4.append(mech_power[3])

# plot RAO wrt wave frequency
plt.plot(T,power_OS1,label='Virginia',color='#CC79A7',marker='*',markersize=14,linewidth=3)
plt.plot(T,power_OS2,label='Frances',color='#009E73',marker='s',markersize=8,linewidth=2,linestyle=(0, (3, 1, 1, 1, 1, 1)))
plt.plot(T,power_PA3,label='Meg',color='#F0E442',marker='>',markersize=12,linewidth=2)
plt.plot(T,power_PA4,label='Jo',color='#56B4E9',marker='o',markersize=12,linewidth=2)

plt.xticks(ticks=T,fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Period [s]',fontsize=20)
plt.ylabel('Power [W]',fontsize=18)
plt.grid()
plt.legend(fontsize=15, markerscale=1,loc='center')
plt.tight_layout()
plt.savefig('array_power.pdf')
plt.clf()

# # plot the PAs
# plt.plot(T,RAO_PA3,label='Meg: Model',color='#F0E442',marker='>',markersize=12,linewidth=2)
# plt.plot(T,RAO_PA4,label='Jo: Model',color='#56B4E9',marker='o',markersize=12,linewidth=2)

# plt.xticks(ticks=T,fontsize=15)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.xlabel('Period [s]',fontsize=20)
# plt.ylabel('RAO [m/m]',fontsize=18)
# plt.grid()
# plt.legend(fontsize=15, markerscale=1)  #,loc='lower left')
# plt.tight_layout()
# plt.savefig('PA_RAOs.pdf')
# plt.clf()

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

# print('Virginia:')
# print('added mass',add_OS1)
# print('damping',damp_OS1)
# print('ex force',ex_OS1)
# print('Frances:')
# print('added mass',add_OS2)
# print('damping',damp_OS2)
# print('ex force',ex_OS2)
# print('Meg:')
# print('added mass',add_PA3)
# print('damping',damp_PA3)
# print('ex force',ex_PA3)
# print('Jo:')
# print('added mass',add_PA4)
# print('damping',damp_PA4)
# print('ex force',ex_PA4)