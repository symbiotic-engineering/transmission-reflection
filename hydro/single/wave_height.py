def wave_height(total, incoming_fse,lam,res):
    import numpy as np
    from scipy.signal import find_peaks
    import matplotlib.pyplot as plt
    import pandas as pd
    import warnings
    warnings.filterwarnings("ignore", category=np.ComplexWarning)
    ##################################################################
    # Extract relevant columns
    z_up = total[int((res/2)*lam-1),:int((res/2)*lam)]
    zinc_up = incoming_fse[int((res/2)*lam-1),:int((res/2)*lam)]
    z_down = total[int((res/2)*lam-1),int((res/2)*lam):]
    zinc_down = incoming_fse[int((res/2)*lam-1),int((res/2)*lam):]

    def calculate_average_wave_height(data):
        peaks, _ = find_peaks(data, height=0)
        troughs, _ = find_peaks(-data, height=0)

        closest_diff = float('inf')
        closest_pair = None
        H_values = []

        for pk_index in peaks:
            for tr_index in troughs:
                diff = abs(data[pk_index] - data[tr_index])
                if diff < closest_diff:
                    closest_diff = diff
                    closest_pair = (pk_index, tr_index)

        if closest_pair:
            peak_index, trough_index = closest_pair
            H = data[peak_index] - data[trough_index]
            H_values.append(H)

        average_H = np.mean(H_values) if H_values else 0
        return average_H

    avg_H_zup = calculate_average_wave_height(z_up)
    avg_H_zincup = calculate_average_wave_height(zinc_up)
    avg_H_zdown = calculate_average_wave_height(z_down)
    avg_H_zincdown = calculate_average_wave_height(zinc_down)

    ref_H = abs(avg_H_zup/avg_H_zincup) - 1
    trans_H = abs(avg_H_zdown/avg_H_zincdown)
    EB1 = trans_H**2 + ref_H**2     # energy balnce
    EB2 = (trans_H + ref_H)**2      # energy balance
    return ref_H, trans_H, EB1, EB2