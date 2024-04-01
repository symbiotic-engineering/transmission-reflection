import sys
import os
# Get the current directory of the script
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
hydro_dir = os.path.join(parent_dir, 'hydro', 'single')
swan_dir = os.path.join(parent_dir, 'swan')

sys.path.append(hydro_dir)
sys.path.append(swan_dir)

import chicken
import alpaca
import numpy as np

w = np.array([1.047])   # wave frequency

Kr_H, Kt_H, Kr_K, Kt_K = chicken.singlebody(w)
KR = [Kr_K[0], Kr_K[0], Kr_K[0], Kr_K[0], Kr_K[0], Kr_K[0]]
KT = [Kt_K[0], Kt_K[0], Kt_K[0], Kt_K[0], Kt_K[0], Kt_K[0]]

alpaca.swanrun(KR,KT)