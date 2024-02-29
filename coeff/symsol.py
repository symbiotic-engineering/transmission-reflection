import numpy as np

def parameters(B, w, M, A, K, e, phi_down, phi_up):
    Z_m = B + 1j*(-w*(M + A) + (K/w))       # step 2: solve eq (11) for complex mechanical impedance
    v = e/Z_m                               # step 3: solve for complex velocity amplitude [eq (10)]
    A_down = phi_down/v                     # step 4: plug v into eq (3) to solve for complex A_+ (downstream amplitude)
    A_up = phi_up/v
    check_sym = A_up - A_down              # step 5: check relationship between complex potential amplitudes
    return A_down, A_up, check_sym, Z_m

def linear(e, A_down, Z_m, A_up):
    import numpy as np
    # solving system of linear equations (symmetrical mode of motion)
    P = (e*A_down)/Z_m                  # define for simplicity
    S = np.array([[-A_down], [P], [P], [0]])
    Q = np.array([[np.conj(A_down), np.conj(A_up), 0, 0],
                  [-1, 0, 0, 1],
                  [0, -1, 1, 0],
                  [-1, 1, -1, 1]])
    invQ = np.linalg.pinv(Q)            # psuedoinverse
    result = np.matmul(invQ, S)
    t, r, R, T = result[:, 0]
    return t, r, R, T

def coefficients(phi_rad, B, w, M, A, K, e, lam):
    t = []
    r = []
    R = []
    T = []
    for i in range (0,2*lam-1):
        phi_down = phi_rad[2*lam-1, 2*lam - i]      # starting at body and moving in -x direction
        phi_up = phi_rad[2*lam-1, 2*lam + i]        # starting at body and moving in +x direction
        A_down, A_up, check_sym, Z_m = parameters(B, w, M, A, K, e, phi_down, phi_up)
        # print('check = zero:', abs(check_sym)) # for symmetric modes of motion this should be zero
        ti, ri, Ri, Ti = linear(e, A_down, Z_m, A_up)  # solve system of equations for coefficients
        t.append(abs(ti))
        r.append(abs(ri))
        R.append(abs(Ri))
        T.append(abs(Ti))
    return t, r, R, T

def find_peaks(data, window_size):
    peaks = []
    for i in range(20, len(data) - 1):
        if data[i] == max(data[max(0, i - window_size):min(len(data), i + window_size + 1)]):
            peaks.append(i)
    return peaks

def avg_peak_height(data,window_size):
    peaks = find_peaks(data,window_size)
    peak_values = data[peaks]
    average_peak_height = np.mean(peak_values)
    return average_peak_height

def avg_coeff(t, r, R, T, lam):
    t = np.array(t)
    r = np.array(r)
    R = np.array(R)
    T = np.array(T)
    
    window_size = int(lam//2.5)
    avg_R = avg_peak_height(R,window_size)
    avg_T = avg_peak_height(T,window_size)
    avg_r = avg_peak_height(r,window_size)
    avg_t = avg_peak_height(t,window_size)
    return avg_R, avg_T, avg_r, avg_t