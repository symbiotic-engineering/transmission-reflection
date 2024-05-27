def postpro(sfgrid_dat,xgrid,ygrid,mxc,myc):
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    import pandas as pd
    import numpy as np

    wecx = [1490 +20, 1590 +20, 1690+20, 1420+20, 1520+20, 1620+20]    # x-position of bodies #, 1790+20,1890+20,1990+20] 
    wecy = [4500, 4500, 4500, 4400, 4400, 4400]    # y-position of row 1 #, 4500, 4500, 4500]
    h_s = pd.read_table(sfgrid_dat, sep="\s+", header=None)
    waveHeight = np.savetxt('/mnt/c/Users/ov162/transmission-reflection/data/atten_1d.csv',h_s,delimiter=',')
    nx, ny = (mxc + 1, myc + 1)
    x = np.linspace(0, xgrid, nx)
    y = np.linspace(0, ygrid, ny)
    X, Y = np.meshgrid(x,y)
    fig,ax = plt.subplots(1,1)
    cp = ax.contourf(X,Y,h_s,40)
    plt.scatter(wecx,wecy,marker='o',color='red')
    ax.set_xlabel('x [m]',fontsize=17)
    ax.set_ylabel('y [m]',fontsize=17)

    ax.ticklabel_format(axis='both',style='sci')
    cbar = fig.colorbar(cp)
    cbar.set_label('Wave Height',fontsize=17,rotation=270,labelpad=17)
    cbar.ax.tick_params(labelsize=15)
    cp.set_clim(0.684, 0.81)
    plt.xticks(fontsize=14, rotation=90)
    plt.yticks(fontsize=14, rotation=90)
    plt.tight_layout()
    plt.savefig('swanfield.pdf')
    plt.show()
    return waveHeight

# import numpy as np
# xgrid = 3000                                # size of grid in x-direction
# ygrid = 5000                                # size of grid in y-direction
# mxc = 300                                   # number of grid points in x (-1)
# myc = 300                                   # number of grid points in y (-1)
# sfgrid_dat = '/mnt/c/Users/ov162/transmission-reflection/run_all/sfgrid.dat'

# waveheight = postpro(sfgrid_dat,xgrid,ygrid,mxc,myc)