def generate_swan_input(KR, KT, d, x, ya, yb):
    import subprocess
    # Calculate xe values
    xe = [xi + d for xi in x]
    # Define the commands
    commands = [
        "PROJ 'southfork' 'A11'",
        "CGRID 0. 0. 0. 3000. 3000. 300 300 CIRCLE 100 0.05 0.25 40",
        "INPGRID BOTTOM 0. 0. 0. 10 10 300. 300.", 
        "READINP BOTTOM -1. 'bathymetry.bot' 1 0 FREE",
        "BOU SHAP JONSWAP 0.77 PEAK DSPR DEGREES",
        "BOU SIDE N CONSTANT PAR 1.3832 6 270 15",
        "BOU SIDE W CONSTANT PAR 1.3832 6 270 15",
        "BOU SIDE E CONSTANT PAR 1.3832 6 270 15",
        "GEN3",
        "OFF WCAP",
        "OFF QUAD",
        "OFF BREA",
    ]

    # Add obstacle lines with varying KT and KR values
    for i in range(len(x)):
        line = f"OBSTACLE TRANS {KT[i]} REFL {KR[i]} RDIFF 1 LINE {x[i]} {ya if i < 3 else yb} {xe[i]} {ya if i < 3 else yb}"
        commands.append(line)

    # Add the remaining commands
    commands.extend([
        "DIFFRAC",
        "FRAME 'SFG' 0 0 0 3000 3000 100 100",
        "OUTPUT OPTIONS '#' BLOCK 4 101",
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