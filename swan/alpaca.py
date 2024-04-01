# script to run swan and post-processing and prepped to talk to capytaine
import numpy as np
import run_swan
import post_process
import run_swan

# Example usage
KR = [0.7, 0.7, 0.7, 0.7, 0.7, 0.7]
KT = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
d = 30
x = [1490, 1590, 1690, 1420, 1520, 1620]
ya = 2000
yb = 900

sfgrid_dat, sfgrid_tbl = run_swan.generate_swan_input(KR, KT, d, x, ya, yb)
post_process.postpro(sfgrid_dat)