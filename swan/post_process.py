def postpro(sfgrid_dat,xgrid,ygrid,mxc,myc):
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    import pandas as pd
    import numpy as np

    #six body reg
    #wecx=[1290+15,1390+15,1490+15,1590+15,1690+15,1790+15]
    #wecy=[4500, 4500, 4500,4500,4500,4500]
    #three body stag: wecx = [1490+15,1590+15,1690+15] wecy=[4450,4500,4450]
    #three-body reg: wecx=[1490+15,1590+15,1690+15] wecy=[4500,4500,4500]
    #six body stag: 
    wecx=[1490+15, 1590+15, 1690+15, 1420+15, 1520+15, 1620+15]
    wecy=[4500,4500,4500,4400,4400,4400]

    h_s = pd.read_table(sfgrid_dat, sep="\s+", header=None)
    #waveHeight = np.savetxt('atten_1d.csv',h_s,delimiter=',')
    waveHeight = np.savetxt('/mnt/c/Users/ov162/transmission-reflection/data/OS_1d.csv',h_s,delimiter=',')
    #h_s = np.loadtxt('/mnt/c/Users/ov162/transmission-reflection/data/PA_3b.csv',delimiter=',')
    nx, ny = (mxc + 1, myc + 1)
    x = np.linspace(0, xgrid, nx)
    y = np.linspace(0, ygrid, ny)
    X, Y = np.meshgrid(x,y)
    fig,ax = plt.subplots(1,1)
    cp = ax.contourf(X,Y,h_s,40)
    plt.scatter(wecx,wecy,marker='_',color='red',s=25, linewidth=2)
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
    print('yip')
    plt.savefig('OS_1d.pdf')
    print('ee')
    plt.show()
    return waveHeight

# import numpy as np
# xgrid = 3000                                # size of grid in x-direction
# ygrid = 5000                                # size of grid in y-direction
# mxc = 300                                   # number of grid points in x (-1)
# myc = 300                                   # number of grid points in y (-1)
# sfgrid_dat = '/mnt/c/Users/ov162/transmission-reflection/run_all/sfgrid.dat'

# waveheight = postpro(sfgrid_dat,xgrid,ygrid,mxc,myc)