This folder will run Capytaine through SWAN start to finish via the big_run.py file.

To run SWAN, you need to have the software installed on your device. Instructions 
are found here: https://gitlab.tudelft.nl/citg/wavemodels/swan

For this script, you will need a bathymetry file, an input file, and initialization file, and 
the SWAN executable. They are all found in this folder as bathymetry.bot, swan_input.swn, 
swaninit, and swan.exe, respectively. The sfgrid.dat and sfgrid.tbl are required for the
post-processing requests in the run_swan function in the swan folder. The Errfile, PRINT file, 
and norm_end file are generated after each run. They can be deleted. (NOTE: the swan_input.swn
file is generated in the swan folder by the run_swan.py script, but it needs to physically
exist in the run_all folder first.)

This folder will first run the near-field Capytaine simulation for the chosen body,
configuration, and sea state using the run_coeffs.wec_run() function. It will compute the 
transmission and reflection coefficients for each body and assign them to the proper obstacle 
in SWAN. The first coefficient in the Kt and Kr arrays corresponds to body #1 and so on 
(as defined in the paper). If a singular body "breaks" energy balance (in the array case, picks 
up dissipated energy from neighboring  bodies), the script will deem that body a non-obstacle.

The run_swan.generate_swan_input() function will generate the input file. The post_process.postpro()
function will generate a .csv file with your data and a figure. The post_pro folder within the
run_all folder manipulates those datasets to generate wave height reduction as a function 
of distance comparison plots. This file is up to your discretion on how to plot and 
present your data.