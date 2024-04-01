def postpro(sfgrid_dat):
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    import pandas as pd
    import numpy as np

    h_s = pd.read_table(sfgrid_dat, sep="\s+", header=None)
    nx, ny = (101, 101)
    x = np.linspace(0, 3000, nx)
    y = np.linspace(0, 3000, ny)
    X, Y = np.meshgrid(x,y)
    fig,ax = plt.subplots(1,1)
    cp = ax.contourf(X,Y,h_s,100)
    ax.set_xlabel('x [m]',fontsize=17)
    ax.set_ylabel('y [m]',fontsize=17)

    ax.ticklabel_format(axis='both',style='sci')
    cbar = fig.colorbar(cp)
    cbar.set_label('Significant Wave Height',fontsize=17,rotation=270,labelpad=17)
    cbar.ax.tick_params(labelsize=15)
    plt.xticks(fontsize=14, rotation=90)
    plt.yticks(fontsize=14, rotation=90)
    plt.tight_layout()
    plt.savefig('test.pdf')
    plt.show()
    return