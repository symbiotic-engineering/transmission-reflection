def generate_swan_input(KR, KT, d, x, ya, yb, H, T, xgrid, ygrid, mxc, myc):
    import subprocess
    # Calculate xe values
    xe = [xi + d for xi in x]
    # Define the commands
    commands = [
        "PROJ 'southfork' 'A11'",                                                   # define project
        # CIRCLE {mdc (meshes in theta)} {flow (lowest discrete frequency[Hz])} {fhigh (highest discrete frequency [Hz])} {msc (grid res in freq space, one less than # of freq)}
        f"CGRID 0. 0. 0. {xgrid} {ygrid} {mxc} {myc} CIRCLE 200 0.06 0.2 20",      # initiate wave propagation grid
        f"INPGRID BOTTOM 0. 0. 0. 10 10 {mxc} {myc}",                               # initiate bottom grid for bathymetry
        "READINP BOTTOM -1. 'bathymetry.bot' 1 0 FREE",                             # define bathymetry (input file from location)
        "WIND 4.92 243",                                                            # define wind input (wind on, wind speed, wind direction)
        #"BOU SHAP JONSWAP 1.54 PEAK DSPR DEGREES",                                  # define JONSWAP spectral shape (alternate value: 0.77)
        "BOU SHAP PM PEAK DSPR DEGREES",                                            # define Pierson-Moskowitz spectrum for deep water waves
        f"BOU SIDE N CONSTANT PAR {H} {T} 270 15",                                  # boundary conditions northern boundary
        f"BOU SIDE S CONSTANT PAR {H} {T} 270 15",                                 # boundary conditions southern boundary
        f"BOU SIDE W CONSTANT PAR {H} {T} 270 15",                                 # boundary conditions western boundary
        f"BOU SIDE E CONSTANT PAR {H} {T} 270 15",                                 # boundary conditions eastern boundary
        "DIFFRAC",                                                                 # energy dissipation due to diffraction ON
        "FRICTION JON CONSTANT",                                                   # energy dissipation due to bottom friction ON
        "PROP BSBT",
        "GEN3 WESTH",                                                               # run generation 3 wave model
        "WCAP",                                                                    # whitecapping energy dissipation ON
        "QUAD",                                                                    # quad wave interactions ON
        "OFF BREA"                                                                  # breaking wave energy dissipation OFF
    ]

    # Add obstacle lines with varying KT and KR values
    for i in range(len(x)):
        # ensuring energy balance
        if KR[i] < 0:
            KR[i] = 0
        EB = (KT[i])**2 + (KR[i])**2
        if KT[i] > 1:
            KT[i] = 1
            KR[i] = 0

        y_start = y_end = [ya,ya,yb,yb,
                           ya,ya,yb,yb,
                           ya,ya,yb,yb,
                           ya,ya,yb,yb,
                           ya,ya,yb,yb]
 
        line = f"OBSTACLE TRANS {KT[i]} REFL {KR[i]} RDIFF 1 LINE {x[i]} {y_start[i]} {xe[i]} {y_end[i]}"
        commands.append(line)

    # Add the remaining commands
    commands.extend([
        f"FRAME 'SFG' 0 0 0 {xgrid} {ygrid} {mxc} {myc}",       # instead of mxc and myc it was 100
        f"OUTPUT OPTIONS '#' BLOCK 4 {mxc + 1}",
        "BLOCK 'SFG' NOHEAD 'sfgrid.dat' LAY 4 HSIGN",
        "BLOCK 'SFG' HEAD 'sfgrid.tbl' HSIGN",
        "TEST 1,0",
        "COMPUTE",
        "STOP",
    ])

    # Write the commands to a file
    with open('swan_input.swn', 'w') as f:
        f.write('\n'.join(commands))

    subprocess.run(['./swan.exe', 'swan_input.swn'])

    # Return the input file name
    return 'sfgrid.dat', 'sfgrid.tbl'