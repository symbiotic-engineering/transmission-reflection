import numpy as np
import matplotlib.pyplot as plt
import bodies
import solve

B = 30*(2*np.pi/180)                                    # wave heading [rad]
depth = 1.37                                            # water depth at basin [m]
T = np.array([1.0,1.2,1.25,1.33,1.39])  #1.003,1.254,1.3])        # wave period from exp [s] (only 1.0,1.25,1.33 for PTO engaged tests. 1.0,1.2,1.25,1.33,1.39 for PTO disengaged)
w = (2*np.pi)/T                                         # wave frequency [rad/s]
g = 9.81
k = w**2/g
reactive = False                         # whether controls are engaged
het = True                              # decide if analyzing hetero or homogeneous array
PA = False                                # IF doing homogeneous array, this determines which architecture
res = 10                                 # grid resolution multiplier

array = bodies.initialize(het,PA)             # create meshed array

RAO_OS1, RAO_OS2, RAO_PA3, RAO_PA4 = [], [], [], []
add_OS1, add_OS2, add_PA3, add_PA4 = [], [], [], []
damp_OS1, damp_OS2, damp_PA3, damp_PA4 = [], [], [], []
ex_OS1, ex_OS2, ex_PA3, ex_PA4 = [], [], [], []
BPTO_OS1, BPTO_OS2, BPTO_PA3, BPTO_PA4 = [], [], [], []
power_OS1, power_OS2, power_PA3, power_PA4 = [], [], [], []
power1, power2, power3, power4 = [], [], [], []
wave_amplitude = [[],[],[]]
emp_radiation = [[],[],[]]
error = [[],[],[]]

## empirically derived damping coeffs from isolated device tests
## just the three periods for PTO engaged tests
# B_diffPA = np.array([2.601353059523+11.0361509650696*1j,0.1571712306061+1.39353317578805*1j,0.4569846084097+4.06849014436796*1j])
# B_diffOS = np.array([5.21423039863538+1.81357171210956*1j,2.14502513133324+0.579566593860094*1j,1.76241880789527+0.457753852086662*1j])
## all five periods tested
B_diffPA = np.array([2.601353059523+11.0361509650696*1j,0.937927576643199+6.28643744791731*1j,0.1571712306061+1.39353317578805*1j,0.4569846084097+4.06849014436796*1j,1.4500030088802+14.7873526149951*1j])
B_diffOS = np.array([5.21423039863538+1.81357171210956*1j,2.46409818130205+0.684169279546558*1j,2.14502513133324+0.579566593860094*1j,1.76241880789527+0.457753852086662*1j,1.48200819928462+0.386856344516752*1j])

# empirically derived damping coeffs from array tests, PTO disengaged
# Bd_vir = np.array([28.660523034874206+126.9728540816521*1j,27.20076324837215+81.03140354379195*1j,24.826928807253346+72.14867163750871*1j])
# Bd_fran = np.array([23.291842420661137+135.57949865296885j,29.47379714707714+78.06298269229453*1j,25.50907937179154+69.8867808578375*1j])
# Bd_meg = np.array([1.2266677151860503-93.80138359019993j,-81.22483790518561-124.79505322193616*1j,-104.14035743647037-121.02988347419092*1j])
# Bd_jo = np.array([-18.564006005286327-145.8236468237688j,-117.2861077516863-127.67332609452676*1j,-132.97385798117375-116.30957769390076*1j])
Bd_vir = np.array([28.660523034874206+126.9728540816521*1j,27.924842028513627+88.61001683278265*1j,27.20076324837215+81.03140354379195*1j,24.826928807253346+72.14867163750871*1j,17.463389920514537+67.42201135724427*1j])
Bd_fran = np.array([23.291842420661137+135.57949865296885j,29.172103691844494+85.38048193975598*1j,29.47379714707714+78.06298269229453*1j,25.50907937179154+69.8867808578375*1j,16.83931066768374+65.94098934998279*1j])
Bd_meg = np.array([1.2266677151860503-93.80138359019993j,-63.27161920377039-130.74541394657783*1j,-81.22483790518561-124.79505322193616*1j,-104.14035743647037-121.02988347419092*1j,-98.82567271173829-126.03085146322015*1j])
Bd_jo = np.array([-18.564006005286327-145.8236468237688j,-82.16213144252684-126.85806877585861*1j,-117.2861077516863-127.67332609452676*1j,-132.97385798117375-116.30957769390076*1j,-108.10513921859607-117.18251379263745*1j])

# ## for PTO disengaged
# megRAOexp = np.array([0.819664657814324,0.916166493356739,0.887306516520389,0.863314218351526,0.992867758219807])
# joRAOexp = np.array([0.477702270828799,0.848824917334879,0.680824985210178,0.713990464053677,0.958269910240572])
# virRAOexp = np.array([323.381890238651,361.268452383110,346.366992147017,322.459657109623,385.006028894270])*(np.pi/180)/k
# franRAOexp = np.array([393.634896229757,384.644580468920,338.447782918374,316.329978665478,396.131175641520])*(np.pi/180)/k

if het:
    # # wave heading = 30deg
    # pos1RAOexp = virRAOexp = np.array([222.769626732413,256.606505334896,336.965180715908])*(np.pi/180)/k
    # pos2RAOexp = franRAOexp = np.array([311.445375226608,354.669863214635,279.004011614259])*(np.pi/180)/k
    # pos3RAOexp = megRAOexp = np.array([0.851771121541293,0.867994446947647,0.886233860061067])
    # pos4RAOexp = joRAOexp = np.array([0.0614839466094648,0.861063004974247,0.647146324540716])
    # wave heading = 0deg
    pos1RAOexp = virRAOexp = np.array([323.381890238651,361.268452383110,346.366992147017,322.459657109623,385.006028894270])*(np.pi/180)/k
    pos2RAOexp = franRAOexp = np.array([393.634896229757,384.644580468920,338.447782918374,316.329978665478,396.131175641520])*(np.pi/180)/k
    pos3RAOexp = megRAOexp = np.array([0.819664657814324,0.916166493356739,0.887306516520389,0.863314218351526,0.992867758219807])
    pos4RAOexp = joRAOexp = np.array([0.477702270828799,0.848824917334879,0.680824985210178,0.713990464053677,0.958269910240572])
    ## RAOS for PTO engaged, het 4b, wave heading = 0 deg
    # pos1RAOexp = virRAOexp = np.array([184.572145311627,195.391623861709,165.997190210489])*(np.pi/180)/k
    # pos2RAOexp = franRAOexp = np.array([131.852630012738,153.262515441623,241.872472048376])*(np.pi/180)/k
    # pos3RAOexp = megRAOexp = np.array([0.773495530727839,0.770132236212850,0.700512165810569])
    # pos4RAOexp = joRAOexp = np.array([0.258407324729989,0.106376997648463,0.0307451184177340])
else:
    if PA:
        # wave heading = 0deg
        pos1RAOexp = laurieRAOexp = np.array([0.0858484279663179,0.104651372235746,0.224911950269082,0.344542057299575,0.394580674608162])
        pos2RAOexp = amyRAOexp = np.array([0.725604161066244,0.934533019923435,1.06503104962314,0.933943156099038,1.05057839656693])
        pos3RAOexp = megRAOexp = np.array([0.780763548356339,0.837774926004850,0.942266698276111,0.837049874715825,0.964715709501276])
        pos4RAOexp = joRAOexp = np.array([0.142736070010829,0.219303581422281,0.596658955054921,0.713137619751200,0.868432685646551])
    else:
        # wave heading = 30deg
        pos1RAOexp = virRAOexp = np.array([244.838050596499,279.178478129899,372.760305464555])*(np.pi/180)/k
        pos2RAOexp = franRAOexp = np.array([383.087742687812,447.702974782584,345.981713099520])*(np.pi/180)/k
        pos3RAOexp = gretaRAOexp = np.array([293.214002674612,291.982847384814,311.308345609112])*(np.pi/180)/k
        pos4RAOexp = eloiseRAOexp = np.array([239.088732269719,363.041856905150,293.509886687945])*(np.pi/180)/k
        # wave heading = 0deg
        # pos1RAOexp = virRAOexp = np.array([309.228753689525,376.784474637044,359.001838403097,352.953996324201,418.158940639490])*(np.pi/180)/k
        # pos2RAOexp = franRAOexp = np.array([364.173892733052,438.289408275332,169.063460623613,157.316360530740,462.835877177627])*(np.pi/180)/k
        # pos3RAOexp = gretaRAOexp = np.array([290.682683775272,375.164478737560,356.294556626036,368.080246646999,389.159325698300])*(np.pi/180)/k
        # pos4RAOexp = eloiseRAOexp = np.array([309.634046596711,363.433614062791,240.697997465533,237.876671971100,386.960953379606])*(np.pi/180)/k

## wave gauge locations from experiments:
#xlocs = np.array([7.093,7.101,8.168,8.173,8.697,8.974,8.967,9.926,13.086,13.072,
#                    13.086,15.08,14.946,16.783,16.78,16.779,19.205,19.2,19.201,21.66])
#ylocs = np.array([-1.185,1.172,-0.456,0.434,-0.008,-0.726,0.719,-0.005,-0.866,-0.047,
#                    0.738,-0.011,3.795,-0.801,0.006,0.778,-2.018,-0.018,1.987,-0.003])
xlocs = np.array([13.086,13.072,13.086,15.08,16.783,16.78,16.779]) - (13.086 + 1.55)
ylocs = np.array([-0.866,-0.047,0.738,-0.011,-0.801,0.006,0.778]) - (-0.1410 - 0.50)

## total wave amplitude from experiments
# wg_number = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
# wg_amp_exp = [[0.0298907659520022,0.0293325490786361,0.0354504770413608,0.0347130310618725,0.0297112053552424,0.0342556822972262,0.0338590620754697,0.0340888702772882,0.0350768531201513,0.0336960205635644,0.0323968750533070,0.0314359072171386,0.0302987506781820,0.0307741198971832,0.0263379692392663,0.0277422035257216,0.0285810659426583,0.0289577544142422,0.0274127509491926,0.0306706231971283],
#                [0.0319177556564394,0.0323068593253320,0.0322781532303995,0.0320792120079917,0.0287044820299621,0.0329302178814983,0.0311944857719814,0.0312402485550472,0.0323303347020503,0.0309218963100902,0.0317841866335105,0.0322292191704004,0.0296665174982208,0.0294324729952453,0.0285451991449859,0.0280614039513836,0.0304093528752781,0.0313865652562018,0.0292021845919162,0.0324239458279942],
#                [0.0341112076039184,0.0329492751549949,0.0324958960778093,0.0327922200171275,0.0335606119510341,0.0334086995562756,0.0320585480407474,0.0325326293450511,0.0317046622601134,0.0319120348567075,0.0316107026751978,0.0318513821422231,0.0334720037655097,0.0306529411814004,0.0309335255946475,0.0289485796752484,0.0294825972813363,0.0304595948273194,0.0273741349622001,0.0317952672279515]]
wg_number = np.array([9,10,11,12,14,15,16])
wg_one = [0.0298907659520022,0.0319177556564394,0.0341112076039184]
wg_amp_exp = [[0.0350768531201513,0.0336960205635644,0.0323968750533070,0.0314359072171386,0.0307741198971832,0.0263379692392663,0.0277422035257216],
                [0.0323303347020503,0.0309218963100902,0.0317841866335105,0.0322292191704004,0.0294324729952453,0.0285451991449859,0.0280614039513836],
                [0.0317046622601134,0.0319120348567075,0.0316107026751978,0.0318513821422231,0.0306529411814004,0.0309335255946475,0.0289485796752484]]


for i in range(np.size(w)):    
    print('running wave period: ',T[i])                                   
    RAO, diff_result, rad_result, ex_force, added_mass, damping, B_PTOemp, mech_power, dissipative_power = solve.hydro(array,B,depth,w[i],reactive,
                                                                            B_diffPA[i],B_diffOS[i],pos3RAOexp[i],
                                                                            pos4RAOexp[i],pos1RAOexp[i],pos2RAOexp[i],
                                                                            Bd_vir[i],Bd_fran[i],Bd_meg[i],Bd_jo[i])   # compute RAOs, excitation force

    ## time to compute the wave elevationnnnnnnnnn
    # total, incoming_fse, grid, radiation, diffraction, elevation_at_gauges = solve.elevation(res,diff_result,rad_result,RAO,xlocs,ylocs,wg_one[i])
    
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

    power1.append(dissipative_power[0])
    power2.append(dissipative_power[1])
    power3.append(dissipative_power[2])
    power4.append(dissipative_power[3])

#     #emp_rad = wg_amp_exp[i] - (np.abs((elevation_at_gauges)))# + (emp_diffraction[i] / wg_diff_exp[i][0])*wg_amp_exp[i][0])
#     #emp_correction = (np.abs((elevation_at_gauges)) + emp_rad) #+ (emp_diffraction[i] / wg_diff_exp[i][0])*wg_amp_exp[i][0])
#     wave_amplitude[i].append(np.abs(elevation_at_gauges))
#     #emp_radiation[i].append(emp_rad)

#     percent_error = -1*((wg_amp_exp[i] - np.abs(elevation_at_gauges))/wg_amp_exp[i])*100
#     error[i].append(percent_error)
#     #print('error',percent_error)

#     plt.figure(figsize=(9, 6))
#     cm = plt.colormaps.get_cmap('seismic')
#     from matplotlib import colors
#     plt.grid()
#     sc = plt.scatter(xlocs, ylocs, c=percent_error,s=400, cmap=cm, 
#                 norm=colors.TwoSlopeNorm(vcenter=0.))

#     x_start = 13.086 + 1.55
#     y_start = -0.1410 - 0.50
#     Dx=Dy=40/50
#     wecx = [[0,0],[0+Dx, 0+Dx]] #[x_start, x_start],[x_start+Dx, x_start+Dx]]
#     wecy = [[0,0+Dy],[0 + (1/2)*Dy,0 + (3/2)*Dy]] #[y_start, y_start + Dy], [y_start + (1/2)*Dy,y_start + (3/2)*Dy]]
#     plt.scatter(wecx[0],wecy[0],marker="|",s = 800, c = 'c')
#     plt.scatter(wecx[1],wecy[1],marker="o",s = 400, c = 'c')

#     for n, txt in enumerate(percent_error):
#         plt.annotate('{0:.2f}'.format(txt), (xlocs[n], ylocs[n]),xytext=(-12, 12),textcoords='offset points',weight="bold",fontsize=14)

    
#     plt.colorbar(sc)
#     plt.xticks(fontsize=15)
#     plt.yticks(fontsize=15)
#     plt.xlabel('x [m]',fontsize=20)
#     plt.ylabel('y [m]',fontsize=20)
#     #plt.ylim([-0.95,0.95])
#     periodnames = ['1.0 s','1.25 s','1.33 s']
#     plt.title('Model Error for ' + periodnames[i],fontsize=20)
#     plt.tight_layout()
#     names = ['1.0map','1.25map','1.33map']
#     plt.savefig(names[i] + '.pdf')
#     plt.clf()

#     # sc = plt.scatter(xlocs, ylocs, c=wg_amp_exp[i], vmin=0.020,vmax=0.050,s=300, cmap=cm)
#     # plt.colorbar(sc)
#     # plt.xticks(fontsize=15)
#     # plt.yticks(fontsize=15)
#     # plt.xlabel('x [m]',fontsize=20)
#     # plt.ylabel('y [m]',fontsize=20)
#     # plt.title('Experiments' + periodnames[i],fontsize=20)
#     # plt.grid()
#     # plt.tight_layout()
#     # names = ['1.0exp','1.25exp','1.33exp']
#     # plt.savefig(names[i] + '.pdf')
#     # plt.clf()

#     # # plot experimental wave amplitude with predicted wave amplitude
#     # plt.plot(wg_number,np.abs(elevation_at_gauges),label='Model',color='#CC79A7',marker='*',markersize=14,linewidth=3)
#     # plt.plot(wg_number,wg_amp_exp[i],label='Experiments',color='#009E73',marker='s',markersize=8,linewidth=2,linestyle=(0, (3, 1, 1, 1, 1, 1)))

#     # plt.xticks(ticks=wg_number,fontsize=15)
#     # plt.xticks(fontsize=15)
#     # plt.yticks(fontsize=15)
#     # plt.xlabel('Wave Gauge',fontsize=20)
#     # plt.ylabel('Wave Amplitude [m]',fontsize=18)
#     # plt.grid()
#     # plt.legend(fontsize=15, markerscale=1)  #,loc='center')
#     # plt.tight_layout()
#     # names = ['1.0wg','1.25wg','1.33wg']
#     # plt.savefig(names[i] + '.pdf')
#     # plt.clf()

#     # plot wave elevation
#     Z = np.abs(np.abs((total)))
#     X = grid[0]
#     Y = grid[1]
#     pcm = plt.pcolormesh(X, Y, Z)
#     plt.xlabel("x")
#     plt.ylabel("y")
#     colorbar = plt.colorbar()
#     #colorbar.set_label(r"Total Wave Elevation, $\eta$")
#     plt.scatter(xlocs,ylocs,marker = 'o', color = 'red', s = 35)
#     colorbar.set_label(r"Wave Amplitude, $\eta$")
#     #pcm.set_clim([0, 2])
#     plt.tight_layout()
#     print('tip')
#     names = ['1.0amp','1.25amp','1.33amp']
#     plt.savefig(names[i] + '.pdf')
#     print('top')
#     plt.clf()

# # plot all experimental results to find trend
# #print('empirical diffraction correction: ',emp_diffraction)

# plt.plot(wg_number,error[0][0],label='T = 1.00 s',color='#CC79A7',marker='*',markersize=14,linewidth=3)
# plt.plot(wg_number,error[1][0],label='T = 1.25 s',color='#009E73',marker='s',markersize=8,linewidth=2,linestyle=(0, (3, 1, 1, 1, 1, 1)))
# plt.plot(wg_number,error[2][0],label='T = 1.33 s',color='#F0E442',marker='>',markersize=12,linewidth=2)

# plt.xticks(ticks=wg_number,fontsize=15)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.xlabel('Wave Gauge',fontsize=20)
# plt.ylabel('Error [%]',fontsize=18)
# plt.grid()
# plt.legend(fontsize=15, markerscale=1)#,loc='center')
# plt.tight_layout()
# plt.savefig('percent_error.pdf')
# # plt.clf()

# plot RAO wrt wave frequency
plt.plot(T,power1,label='Virginia',color='#CC79A7',marker='*',markersize=14,linewidth=3)
plt.plot(T,power2,label='Frances',color='#009E73',marker='s',markersize=8,linewidth=2,linestyle=(0, (3, 1, 1, 1, 1, 1)))
plt.plot(T,power3,label='Meg',color='#F0E442',marker='>',markersize=12,linewidth=2)
plt.plot(T,power4,label='Jo',color='#56B4E9',marker='o',markersize=12,linewidth=2)

plt.xticks(ticks=T,fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Period [s]',fontsize=20)
plt.ylabel('Radiated Power [W]',fontsize=18)
plt.grid()
plt.legend(fontsize=15, markerscale=1)
plt.tight_layout()
plt.savefig('test.pdf')
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

# note TO OLIVIA: (np.abs(total)/np.abs(incoming_fse)) - np.abs(total/incoming_fse) = numerically zero
# note TO OLIVIA: when modeling scaled up vs scaled down verisons, the PAs produce the same RAO, 
# but the OSWEC RAO changes VERY SIGNIFICANTLY, from around 4 to aroun 0.08