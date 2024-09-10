def generate_swan_input(KR, KT, d, x, ya, yb, H, T, xgrid, ygrid, mxc, myc, attenuator, six, swan_stag):
    import subprocess
    # Calculate xe values
    xe = [xi + d for xi in x]
    # Define the commands
    commands = [
        "PROJ 'southfork' 'A11'",
        f"CGRID 0. 0. 0. {xgrid} {ygrid} {mxc} {myc} CIRCLE 100 0.05 0.25 40",
        f"INPGRID BOTTOM 0. 0. 0. 10 10 {mxc} {myc}", 
        "READINP BOTTOM -1. 'bathymetry.bot' 1 0 FREE",
        "WIND 4.92 243",
        "BOU SHAP JONSWAP 1.54 PEAK DSPR DEGREES",       #0.77
        f"BOU SIDE N CONSTANT PAR {H} {T} 270 15",
        f"BOU SIDE S CONSTANT PAR {H} {T} 270 15",
        f"BOU SIDE W CONSTANT PAR {H} {T} 270 15",
        f"BOU SIDE E CONSTANT PAR {H} {T} 270 15",
        "DIFFRAC",
        "FRICTION JON CONSTANT",
        "PROP BSBT",
        "GEN3 WESTH",
        "WCAP",
        "QUAD",
        "OFF BREA",
    ]

    # Add obstacle lines with varying KT and KR values
    for i in range(len(x)):
        # ensuring energy balance
        if KR[i] < 0:
            KR[i] = 0
        EB = (KT[i])**2 + (KR[i])**2
        if EB > 1:
            KT[i] = 1
            KR[i] = 0
        # indexing due to how problem was formulating
        if six:
            use_condition = i < 3
        else:
            use_condition = i == 1 if swan_stag else i < 3
        # attenuator has more length than other devices
        if attenuator:
            length = 29
            y_start = ya if use_condition else yb
            y_end = ya + length if use_condition else yb - length
        else:
            y_start = ya if use_condition else yb
            y_end = ya if use_condition else yb
 
        line = f"OBSTACLE TRANS {KT[i]} REFL {KR[i]} RDIFF 1 LINE {x[i]} {y_start} {xe[i]} {y_end}"
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