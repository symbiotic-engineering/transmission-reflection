import numpy as np
import jonswap
import matplotlib.pyplot as plt
from scipy.integrate import simps

# for single OSWEC, reactive controls
w = np.array([0.7,0.8,0.9,1.0,1.1,1.25,1.3])
kt = np.array([0.8111207754147617,0.7733889383923015,0.7707245622136073,0.7820227405424823,0.7623678809417559,0.6391950252758175,0.5852990549054277])
kr = np.array([0.2829520828664469,0.22748709933258504,0.2044217536608628,0.21996551735827952,0.2585265888662087,0.3426032118086353,0.3542309382683648])

S_j = jonswap.get_spectra(w)
max_Sj = np.max(S_j)
norm_Sj = S_j / max_Sj

kt_pdf = []
kr_pdf = []

for i in range(np.size(kt)):
    k_t_pdf = norm_Sj[i] * kt[i]
    k_r_pdf = norm_Sj[i] * kr[i]

    kt_pdf.append(k_t_pdf)
    kr_pdf.append(k_r_pdf)

expected_kt = np.sqrt(simps(kt,kt_pdf))
expected_kr = np.sqrt(simps(kr,kr_pdf))
print('expected Kt = ',expected_kt)
print('expected Kr = ',expected_kr)

plt.plot(w,kt_pdf,label='$K_t$')
plt.plot(w,kr_pdf,label='$K_r$')
plt.plot(w,norm_Sj,label='Normalized JONSWAP')
plt.legend()
plt.xlabel('$\omega$ [rad/s]')
plt.ylabel('Normalized Value')
plt.tight_layout()
plt.savefig('norm_coeffs.pdf')