def disturbance(kd,lam,res):
    import pandas as pd
    import numpy as np

    xR = np.linspace(-1*int(res*lam),0,int((res/2)*lam))
    xD = np.linspace(0,int(res*lam),int((res/2)*lam))

    #########################################################################################################

    # Extract relevant columns
    kd_up = kd[int((res/2)*lam)-1,:int((res/2)*lam)]

    def find_peaks(data, lambda_wave, xR):
        peaks = []
        for i in range(len(data)):
            if i < 0 or i >= len(data):
                continue
            if i == 0 or i == len(data) - 1:
                continue
            if xR[i] < -1*int(lambda_wave) and data[i] > data[i - 1] and data[i] > data[i + 1]:
                peaks.append(data[i])
        return peaks
    def avg_peak_height(data, lambda_wave, xR):
        peak_values = find_peaks(data, lambda_wave, xR)
        average_peak_height = np.mean(peak_values)
        return average_peak_height


    avg_kd_up = avg_peak_height(kd_up, lam, xR)
    ref_K = abs(avg_kd_up) - 1

    ######################################## FOR DOWNSTREAM ###########################################
    kd_down = kd[int((res/2)*lam)-1,int((res/2)*lam):]

    def peaks_down(series, lambda_wave):
        peaks = []
        for i in range(1, len(series) - 1):  
            if np.angle(series[i] - series[i - 1]) > 0 and np.angle(series[i] - series[i + 1]) > 0:
                if xD[i] < 2*lambda_wave:  # Limiting to x < lambda_wave
                    peaks.append(series[i].real)  
        return peaks

    def avg_peak_down(data, lambda_wave):
        peak_values = peaks_down(data, lambda_wave)
        avg_peak_values = np.mean(peak_values)
        return avg_peak_values

    avg_distance_kd_down = avg_peak_down(kd_down, lam)
    trans_K= abs(avg_distance_kd_down)

    EB1 = trans_K**2 + (ref_K - 1)**2       # energy balance check
    EB2 = (trans_K + (ref_K - 1))**2       # energy balance check
    return ref_K, trans_K, EB1, EB2