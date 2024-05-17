def generate_swan_input(KR, KT, d, x, ya, yb, H, T, xgrid, ygrid, mxc, myc):
    import subprocess
    # Calculate xe values
    xe = [xi + d for xi in x]
    # Define the commands
    commands = [
        "PROJ 'southfork' 'A11'",
        f"CGRID 0. 0. 0. {xgrid} {ygrid} {mxc} {myc} CIRCLE 100 0.05 0.25 40",
        f"INPGRID BOTTOM 0. 0. 0. 10 10 {mxc} {myc}", 
        "READINP BOTTOM -1. 'bathymetry.bot' 1 0 FREE",
        "BOU SHAP JONSWAP 0.77 PEAK DSPR DEGREES",
        f"BOU SIDE N CONSTANT PAR {H} {T} 270 15",
        f"BOU SIDE S CONSTANT PAR {H} {T} 270 15",
        f"BOU SIDE W CONSTANT PAR {H} {T} 270 15",
        f"BOU SIDE E CONSTANT PAR {H} {T} 270 15",
        "GEN3",
        "OFF WCAP",
        "OFF QUAD",
        "OFF BREA",
    ]

    # Add obstacle lines with varying KT and KR values
    for i in range(len(x)):
        if KT[i] > 1:
            KT[i] = 1 - KR[i]
        line = f"OBSTACLE TRANS {KT[i]} REFL {KR[i]} RDIFF 1 LINE {x[i]} {ya if i < 3 else yb} {xe[i]} {ya if i < 3 else yb}"
        commands.append(line)

    # Add the remaining commands
    commands.extend([
        "DIFFRAC",
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