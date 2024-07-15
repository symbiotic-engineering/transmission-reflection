To use this repository, it is recommended to create an environment using the "requirements.txt" file provided. If not, you will need to install Capytaine (https://capytaine.github.io/stable/) and potentially other python packages. Regardless, you will need to install SWAN software (https://gitlab.tudelft.nl/citg/wavemodels/swan) if you want to run the SWAN files (in "swan" folder and "run_all" folder).

This code finds the transmission (Kt) and reflection coefficients (Kr) for different wave energy converter archetypes and computes wave height reduction as a function of distance from the device. The devices modeled in this code are:

1. Heaving Point Absorber Wave Energy Converter
2. Oscillating Surge Wave Energy Converter
3. Attenuator Wave Energy Converter
4. Floating Breakwater

The transmission and reflection coefficients are calculated in the "hydro" folder using the "run_coeffs.py" script. More details are found in the folder README.md file. The SWAN command file is generated in the "swan" folder. To alter the input commands, see the SWAN documentation (linked above) for details on their functions. The using the "big_run.py" script in the "run_all" folder will do all of the following:

1. Generate the floating body of your choice
2. Solve the hydrodynamics of the floating body and the surrounding wave elevation
3. Compute the transmission and reflection coefficients for each body
4. Input those coefficients as obstacles in a SWAN simulation
5. Compute the wave height on the computational grid due to said obstacles
6. Save the wave height data as a .csv file and generate a contour plot
