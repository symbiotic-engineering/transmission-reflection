# working with the potentials from breakwater limiting case
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# initializing parameters
k = 0.1117            #0.0529          # finite depth, k = 0.1117 in infinite depth
w = 1.047
T = (2*np.pi)/w       # period [s]
A_I = 1               # incoming wave amplitude
g = 9.81
z = 0
xR = np.linspace(-100,0,75)
xD = np.linspace(0,100,75)
xt = np.linspace(-100,100,150)
t = 0
rho = 1025          # fluid density [kg/m^3]
y0 = np.zeros_like(xt)
# all the exponentials
ekz = np.exp(k*z)
eikx_R = np.exp(1j*k*xR)
eikx_D = np.exp(-1*1j*k*xD)
eiwt = np.exp(1j*w*t)

#########################################################################################################

# compare wave elevation from capytaine to calcualted from my potentials
# working with data from breakwater limiting case
capy_elev = "dif_plus_inc_data_bw.csv"
el = pd.read_csv(capy_elev, header=None)
zeta = el.applymap(lambda x: complex(x.replace(' ', '')))
z_up = -1*zeta.iloc[:75,74]
z_down = -1*zeta.iloc[75:,74]

# finding amplitudes from capytaine wave elevation
# diffracted elevation
dif_el = "dif_el_data_bw.csv"
d_el = pd.read_csv(dif_el, header=None, dtype=str)
zeta_d = d_el.applymap(lambda x: complex(x.replace('j', 'j').replace('+-', '-')))
zd_up = -1*zeta_d.iloc[:75,74]
zd_down = -1*zeta_d.iloc[75:,74]
Ad_ref = zd_up*np.exp(-1j*k*xR - 1j*w*t)
Ad_tran = zd_down*np.exp(1j*k*xD - 1j*w*t)

# incident elevation
inc_el = "incoming_fse_data_bw.csv"
inc_el = pd.read_csv(inc_el, header=None)
inc = inc_el.applymap(lambda x: complex(x.replace(' ', '')))
zinc_up = -1*inc.iloc[:75,74]
zinc_down = -1*inc.iloc[75:,74]
Ainc_up = zinc_up*np.exp(1j*k*xR - 1j*w*t)
Ainc_down = zinc_down*np.exp(1j*k*xD - 1j*w*t)

# total
A_tot_up = Ad_ref + Ainc_up
A_tot_down = Ad_tran + Ainc_down

# wave elevation from capytaine
plt.plot(xR,zd_up,color='red',linestyle=':',label='Diffracted Upstream')
plt.plot(xD,zd_down,color='black',linestyle=':',label='Diffracted Downstream')
plt.plot(xR,z_up,color='blue',label='Total')
plt.plot(xD,z_down,color='blue')
plt.plot(xt,y0,color='green',linestyle='--',label='Still Water Line')
plt.legend()
plt.xlabel('x [m]')
plt.ylabel('Wave Elevation [m]')
plt.savefig(f'plots/capy_elev.pdf')
plt.show()

# incident potential based on wave elevation
phi_i_up = ((1j*g*Ainc_up)/w)*ekz*np.exp(-1*1j*k*xR)*eiwt
phi_i_down = ((1j*g*Ainc_down)/w)*ekz*np.exp(-1*1j*k*xD)*eiwt

# diffracted potential data
csv_file_path = "potential_data_bw.csv"
df = pd.read_csv(csv_file_path, header=None)
phi_complex = df.applymap(lambda x: complex(x.replace(' ', '')))
pR_comp = phi_complex.iloc[74,:75]
pD_comp = phi_complex.iloc[74,75:]

# incident potential data
airys = "airy_phi_inc.csv"
airy = pd.read_csv(airys, header=None)
airy_phi = airy.applymap(lambda x: complex(x.replace(' ', '')))
pinc_up = airy_phi.iloc[74,:75]
pinc_down = airy_phi.iloc[74,75:]

# calculating complex amplitude from potential
A_inc_up = (w*pinc_up)/(1j*g*ekz*np.exp(-1*1j*k*xR)*eiwt)
A_inc_down = (w*pinc_down)/(1j*g*ekz*np.exp(-1*1j*k*xD)*eiwt)
A_R = ((w*pR_comp)/(1j*g*ekz*eikx_R*eiwt))
A_D = (w*pD_comp)/(1j*g*ekz*eikx_D*eiwt)
Auptot = A_R + A_inc_up
Adowntot = A_D + A_inc_down

# wave elevation from amplitudes found via potential
E_inc_up = np.real(A_inc_up*np.exp(-1j*k*xR + 1j*w*t))
E_inc_down = np.real(A_inc_down*np.exp(-1j*k*xD + 1j*w*t))
E_R = np.real(A_R*np.exp(1j*k*xR + 1j*w*t))    
E_D = np.real(A_D*np.exp(-1j*k*xD + 1j*w*t))
E_up_tot = E_inc_up + E_R
E_down_tot = E_inc_down + E_D

#########################################################################################################
# time averaged energy flux using amplitudes
mult = (1/4)*rho*(g**2/w)         # multiplier for both up and downstream energy flux
dEdt_up = mult*(np.abs(A_inc_up)**2 - np.abs(A_R)**2)
dEdt_down = mult*(np.abs(A_D)**2)
dEdtinc_up = mult*np.abs(A_inc_up)**2
dEdtinc_down = mult*np.abs(A_inc_down)**2

# plt.plot(xR,dEdt_up,color='black',label='Upstream')
# plt.plot(xD,dEdt_down,color='green',label='Downstream')
# plt.plot(xR,dEdtinc_up,color='black',linestyle=':',label='Incident')
# plt.plot(xD,dEdtinc_down,color='green',linestyle=':')
# plt.legend()
# plt.xlabel('x [m]')
# plt.ylabel('Energy Flux [kg-m/s^3]')
# plt.tight_layout()
# plt.savefig(f'plots/approx_eflux.pdf')
# plt.show()

#########################################################################################################

# time averaged energy flux using capytaine potentials, derivatives, and integration (maha method)
pup_tot = pinc_up + pR_comp
pdown_tot = pinc_down + pD_comp
# plt.plot(xR,pup_tot,color='red',label='pup')
# plt.plot(xD,pdown_tot,color='black',label='pdown')
# plt.xlabel('x [m]')
# plt.ylabel('Energy Flux [kg-m/s^3]')
# plt.tight_layout()
# plt.show()

def z_function(zd):
    # finding the real part of the potential derivatives by dx and dt
    z_func = np.exp(zd*k)

    comp_dphi_R_dx = -1j*k*pup_tot*z_func
    comp_dphi_D_dx = -1j*k*pdown_tot*z_func
    comp_dphi_R_dt = 1j*w*pup_tot*z_func
    comp_dphi_D_dt = 1j*w*pdown_tot*z_func

    dphi_R_dx = np.real(comp_dphi_R_dx)
    dphi_D_dx = np.real(comp_dphi_D_dx)
    dphi_R_dt = np.real(comp_dphi_R_dt)
    dphi_D_dt = np.real(comp_dphi_D_dt)

    # pre-integration (inside the integral over z)
    ins_R = dphi_R_dx*dphi_R_dt   # reflection
    ins_D = dphi_D_dx*dphi_D_dt   # diffraction
    return [ins_R,ins_D]


int_R = np.zeros(75)  # Assuming the size of the array is 75
int_D = np.zeros(75)
# Loop from z = -50 to z = 0
for z_value in range(-50,1):
    # integrating the derivatives over z
    ins_R, ins_D = z_function(z_value)
    int_R += ins_R       # integral of multiplied derivative potentials over z
    int_D += ins_D

# now to find the time average of these values
# intR = int_R/np.exp(2*1j*w*t)
# uptest = (np.exp(2*1j*w*T) - np.exp(2*1j*w*t))/(2*1j*w)
# NR_t = np.real(intR*uptest)
NR_t = np.real((int_R*np.exp(2*1j*w*T)/np.exp(2*1j*w*t)) - int_R/np.exp(2*1j*w*t))
ND_t = np.real(int_D*np.exp(2*1j*w*T)/np.exp(2*1j*w*t))

# time averaged energy flux!
dE_R_dt = (NR_t*-rho)/T
dE_D_dt = (ND_t*-rho)/T

plt.plot(xR,dE_R_dt,color='black',label='Upstream')
plt.plot(xD,dE_D_dt,color='green',label='Downstream')
plt.legend()
plt.xlabel('x [m]')
plt.ylabel('Energy Flux [kg-m/s^3]')
plt.tight_layout()
plt.savefig(f'plots/int_eflux.pdf')
plt.show()

##################################################################################################################
# # compare wave elevations
# plt.plot(xR,E_R,color='red',linestyle=':',label='Dif Up Pot')
# plt.plot(xD,E_D,color='black',linestyle=':',label='Dif Down Pot')
# plt.plot(xR,E_up_tot,color='blue',label='Total Pot')
# plt.plot(xD,E_down_tot,color='blue')
# plt.plot(xR,E_R,color='red',linestyle='-.',label='Dif Up Elev')
# plt.plot(xD,E_D,color='black',linestyle='-.',label='Dif Down Elev')
# plt.plot(xR,E_up_tot,color='blue',linestyle=':',label='Total')
# plt.plot(xD,E_down_tot,color='blue',linestyle=':')
# plt.plot(xR,E_inc_up,color='pink',linestyle=':',label='Incident')
# plt.plot(xD,E_inc_down,color='pink',linestyle=':')
# plt.plot(xt,y0,color='green',linestyle='--',label='Still Water Line')
# plt.legend()
# plt.xlabel('x [m]')
# plt.ylabel('Wave Elevation')
# plt.savefig(f'plots/comp_elev.pdf')
# plt.show()

# # plotting wave elevation based on potentials
# plt.plot(xR,E_R,color='red',linestyle=':',label='Diffracted Upstream')
# plt.plot(xD,E_D,color='black',linestyle=':',label='Diffracted Downstream')
# plt.plot(xR,E_up_tot,color='blue',label='Total')
# plt.plot(xD,E_down_tot,color='blue')
# plt.plot(xt,y0,color='green',linestyle='--',label='Still Water Line')
# plt.legend()
# plt.xlabel('x [m]')
# plt.ylabel('Wave Elevation')
# plt.savefig(f'plots/pot_based_elev.pdf')
# plt.show()

# # comparing capytaine airys wave potential and potential found from elevation
# plt.plot(xR,phi_i_up,color='blue',label='Elev')
# plt.plot(xD,phi_i_down,color='blue')
# plt.plot(xR,pinc_up,color='red',linestyle=':',label='Airys')
# plt.plot(xD,pinc_down,color='red',linestyle=':')
# plt.plot(xt,y0,color='green',linestyle='--',label='Still Water Line')
# plt.legend()
# plt.xlabel('x [m]')
# plt.ylabel('Incident Potential')
# plt.savefig(f'plots/compare_pinc.pdf')
# plt.show()

# # plotting amplitudes from capytaine wave elevation
# plt.plot(xR,np.abs(Ad_ref),color='red',linestyle='-.',label='Reflected')
# plt.plot(xD,np.abs(Ad_tran),color='black',linestyle='-.',label='Transmitted')
# plt.plot(xR,np.abs(Ainc_up),color='blue',linestyle=':',label='Incident')
# plt.plot(xD,np.abs(Ainc_down),color='blue',linestyle=':')
# plt.plot(xR,np.abs(A_tot_up),color='blue',label='Total')
# plt.plot(xD,np.abs(A_tot_down),color='blue')
# plt.plot(xt,y0,color='green',linestyle='--',label='Still Water Line')
# plt.legend()
# plt.xlabel('x [m]')
# plt.ylabel('Amplitude [m]')
# plt.savefig(f'plots/dif_fs_amplitudes.pdf')
# plt.show()

# # plotting capytaine potentials
# plt.plot(xR,pinc_up,color='blue',label='Incident')
# plt.plot(xD,pinc_down,color='blue')
# plt.plot(xR,np.real(pR_comp),color='red',label='Diffracted Upstream')
# plt.plot(xD,np.real(pD_comp),color='black',label='Diffracted Downstream')
# plt.xlabel('x [m]')
# plt.ylabel('Potential [m^2/s]')
# plt.legend()
# plt.tight_layout()
# plt.savefig(f'plots/potential.pdf')
# plt.show()

# # compare diffracted amplitude from potential and from elevation

# plt.plot(xR,np.abs(A_R),color='red',label='Reflected Pot')
# plt.plot(xD,np.abs(A_D),color='black',label='Transmitted Pot')
# plt.plot(xR,np.abs(Ad_ref),color='red',linestyle='-.',label='Reflected Elev')
# plt.plot(xD,np.abs(Ad_tran),color='black',linestyle='-.',label='Transmitted Elev')
# plt.plot(xR,np.abs(A_tot_up),color='blue',linestyle='-.',label='Total')
# plt.plot(xD,np.abs(A_tot_down),color='blue',linestyle='-.')
# plt.plot(xR,np.abs(A_inc_up),color='pink',linestyle='-.',label='Incident')
# plt.plot(xD,np.abs(A_inc_down),color='pink',linestyle='-.')
# plt.plot(xt,y0,color='green',linestyle='--',label='Still Water Line')
# plt.legend()
# plt.xlabel('x [m]')
# plt.ylabel('Amplitude [m]')
# plt.tight_layout()
# plt.savefig(f'plots/compare_amp.pdf')
# plt.show()

# # compare total amplitude from potential and elevation
# plt.plot(xR,A_tot_up,color='blue',label='Total Elev')
# plt.plot(xD,A_tot_down,color='blue')
# plt.plot(xR,Auptot,color='black',linestyle=':',label='Total Pot')
# plt.plot(xD,Adowntot,color='black',linestyle=':')
# plt.plot(xt,y0,color='green',linestyle='--',label='Still Water Line')
# plt.legend()
# plt.xlabel('x [m]')
# plt.ylabel('Amplitude [m]')
# plt.savefig(f'plots/tot_amp_compare.pdf')
# plt.show()

# # plotting amplitude from potentials

# plt.plot(xR,np.real(A_R),color='red',label='Reflected')
# plt.plot(xD,np.real(A_D),color='black',label='Transmitted')
# plt.plot(xt,np.real(A_inc),color='blue',label='Incident')
# plt.plot(xt,y0,color='green',linestyle='--',label='Still Water Line')
# plt.legend()
# plt.xlabel('x [m]')
# plt.ylabel('Amplitude [m]')
# plt.tight_layout()
# plt.savefig(f'plots/amplitude.pdf')
# plt.show()

###############################################################################################################

# # Finding the wave height --> reflection
# xR_lower_thresh = -100
# xR_upper_thresh = -25

# # Find the indices where xR is between the specified thresholds
# indices_within_thresh = np.where((xR >= xR_lower_thresh) & (xR <= xR_upper_thresh))[0]

# # Check if there are any points within the specified range
# if len(indices_within_thresh) > 0:
#     # Find the index of the maximum and minimum values within the specified range
#     max_index_ER = indices_within_thresh[np.argmax(E_R[indices_within_thresh])]
#     min_index_ER = indices_within_thresh[np.argmin(E_R[indices_within_thresh])]

#     # Calculate the maximum and minimum values
#     max_ER = E_R[max_index_ER]
#     min_ER = E_R[min_index_ER]

#     # Calculate the wave height and reflection coefficient
#     H_R = max_ER - min_ER
#     KR = (H_R / (2 * A_I)) - 1

#     print('Reflected wave height:', H_R)
#     print('Reflection Coefficient:', KR)
# else:
#     print("No points within the specified range.")

# # transmission
# xD_thresh = 20
# D_index = np.where(xD >= xD_thresh)[0][0]
# max_index_ED = D_index + np.argmax(E_D[D_index:])
# max_ED = E_D[max_index_ED]
# min_index_ED = D_index + np.argmin(E_D[D_index:])
# min_ED = E_D[min_index_ED]
# H_D = max_ED - min_ED
# KT = H_D/(2*A_I)
# print('Transmitted wave height',H_D)
# print('Transmission Coefficient',KT)

# #plt.plot(xt,E_inc,color='blue',linestyle=':',label='Incident')
# plt.plot(xR,E_R,color='red',label='Upstream')
# plt.plot(xD,E_D,color='black',label='Transmitted')
# plt.plot(xt,y0,color='green',linestyle='--',label='Still Water Line')
# plt.plot(xR[max_index_ER], max_ER, 'go',markersize=6) 
# plt.plot(xR[min_index_ER], min_ER, 'go',markersize=6)  
# plt.plot(xD[max_index_ED], max_ED, 'go',markersize=6)  
# plt.plot(xD[min_index_ED], min_ED, 'go',markersize=6)  
# plt.legend()
# plt.xlabel('x [m]')
# plt.ylabel('Wave Elevation [m]')
# plt.tight_layout()
# plt.savefig(f'plots/elevation.pdf')
# plt.show()

###########################################################################################################

# # Finding the wave height --> reflection
# xr_lt = -100
# xr_ut = -5

# # Find the indices where xR is between the specified thresholds
# indices = np.where((xR >= xr_lt) & (xR <= xr_ut))[0]

# # Check if there are any points within the specified range
# if len(indices) > 0:
#     # Find the index of the maximum and minimum values within the specified range
#     max_index_zu = indices[np.argmax(z_up[indices])]
#     min_index_zu = indices[np.argmin(z_up[indices])]

#     # Calculate the maximum and minimum values
#     max_zu = z_up[max_index_zu]
#     min_zu = z_up[min_index_zu]

#     # Calculate the wave height and reflection coefficient
#     hu = max_zu - min_zu
#     kr = (hu / (2 * A_I)) - 1

#     print('Reflected wave height:', hu)
#     print('Reflection Coefficient:', kr)
# else:
#     print("No points within the specified range.")

# # transmission
# xd_thresh = 5
# zd_index = np.where(xD >= xd_thresh)[0][0]
# max_index_zd = zd_index + np.argmax(z_down[zd_index:])
# print('zd_index:', zd_index)
# print('max_index_zd:', max_index_zd)
# print('Length of z_down DataFrame:', len(z_down))

# max_zd = z_down[max_index_zd]
# min_index_zd = zd_index + np.argmin(z_down[zd_index:])
# min_zd = z_down[min_index_zd]
# hd = max_zd - min_zd
# kt = hd/(2*A_I)
# print('Transmitted wave height',hd)
# print('Transmission Coefficient',kt)

# plt.plot(xR,z_up,color='red',label='Upstream')
# plt.plot(xD,z_down,color='black',label='Downstream')
# #plt.plot(xt,E_inc,color='blue',linestyle=':',label='Incident')
# plt.plot(xt,y0,color='green',linestyle='--',label='Still Water Line')
# # plt.plot(xR[max_index_zu], max_zu, 'go',markersize=6) 
# # plt.plot(xR[min_index_zu], min_zu, 'go',markersize=6)  
# # plt.plot(xD[max_index_zd], max_zd, 'go',markersize=6)  
# # plt.plot(xD[min_index_zd], min_zd, 'go',markersize=6) 
# plt.legend()
# plt.xlabel('x [m]')
# plt.ylabel('Wave Elevation [m]')
# plt.savefig(f'plots/capy_elev.pdf')
# plt.show()
