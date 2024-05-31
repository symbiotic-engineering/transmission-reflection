def swanrun(KR,KT,x,ya,yb,H,T,xgrid,ygrid,mxc,myc,attenuator):
    # script to run swan and post-processing and prepped to talk to capytaine
    import numpy as np
    import run_swan
    import post_process
    import run_swan

    d = 30

    sfgrid_dat, sfgrid_tbl = run_swan.generate_swan_input(KR, KT, d, x, ya, yb, H, T, xgrid, ygrid, mxc, myc,attenuator)
    waveHeight = post_process.postpro(sfgrid_dat,xgrid,ygrid,mxc,myc)
    return waveHeight