def wave_height(total, incoming_fse,lam,xtrans,ytrans,farm,rel_dim,nx,ny,x1,x2,y1,y2):
    import numpy as np
    from scipy.signal import find_peaks
    import pandas as pd
    import matplotlib.pyplot as plt
    import warnings
    warnings.filterwarnings("ignore", category=np.ComplexWarning)
    ##################################################################
    # Extract relevant columns
    zinc_up = incoming_fse[int(ny/2),:int((nx/2)-rel_dim)]
    zinc_down = incoming_fse[int(ny/2),int((nx/2)+rel_dim):]
    z_up = total[int(ny/2),:int((nx/2)-rel_dim)]
    z_down = total[int(ny/2),int((nx/2)+rel_dim):]

    avg_H_zup = np.array([np.mean(abs(z_up))])
    avg_H_zincup = np.array([np.mean(abs(zinc_up))])
    avg_H_zdown = np.array([np.mean(abs(z_down))])
    avg_H_zincdown = np.array([np.mean(abs(zinc_down))])


    if farm == True:
        z_upWEC1 = total[int((ny/2) + (ytrans[0]/2)),:int((nx/2)-rel_dim+(xtrans[0]/2))]
        z_upWEC3 = total[int((ny/2) + (ytrans[1]/2)),:int((nx/2)-rel_dim+(xtrans[1]/2))]

        z_downWEC1 = total[int((ny/2) + (ytrans[0]/2)),int((nx/2)+rel_dim+(xtrans[0]/2)):]
        z_downWEC3 = total[int((ny/2) + (ytrans[1]/2)),int((nx/2)+rel_dim+(xtrans[1]/2)):]

        avg_H_zup = np.array([np.mean(abs(z_upWEC1)),np.mean(abs(z_up)),np.mean(abs(z_upWEC3))])
        avg_H_zdown = np.array([np.mean(abs(z_downWEC1)),np.mean(abs(z_down)),np.mean(abs(z_downWEC3))])

    # # Create the grid
    # grid_x, grid_y = np.meshgrid(np.linspace(x1, x2, nx), np.linspace(y1, y2, ny))
    # plt.figure(figsize=(10, 10))
    # plt.contourf(grid_x, grid_y, total, cmap='viridis', levels=50)
    # plt.colorbar(label='Total Values')

    # # Plot z_up and z_down lines
    # y_level = grid_y[int((ny/2)*lam-1), 0]
    # plt.plot(grid_x[0, :int((nx/2)*lam - rel_dim)], z_up + y_level, 'r-', label='z_up')
    # plt.plot(grid_x[0, int((nx/2)*lam + rel_dim):], z_down + y_level, 'b-', label='z_down')

    # # Plot farm-related lines if farm is True
    # if farm:
    #     y_level_wec1 = grid_y[int((ny/2)*lam-1 + ytrans[0]/2), 0]
    #     y_level_wec3 = grid_y[int((ny/2)*lam-1 + ytrans[1]/2), 0]
    #     plt.plot(grid_x[0, :int((nx/2)*lam + (xtrans[0]/2) - rel_dim)], z_upWEC1 + y_level_wec1, 'g-', label='z_upWEC1')
    #     plt.plot(grid_x[0, :int((nx/2)*lam + (xtrans[1]/2) - rel_dim)], z_upWEC3 + y_level_wec3, 'm-', label='z_upWEC3')

    #     plt.plot(grid_x[0, int((nx/2)*lam + (xtrans[0]/2) + rel_dim):], z_downWEC1 + y_level_wec1, 'y-', label='z_downWEC1')
    #     plt.plot(grid_x[0, int((nx/2)*lam + (xtrans[1]/2) + rel_dim):], z_downWEC3 + y_level_wec3, 'c-', label='z_downWEC3')

    # plt.xlabel('X coordinates')
    # plt.ylabel('Y coordinates')
    # plt.legend()
    # plt.title('Geographical Plot of z_up and z_down Lines')
    # plt.show()


    ref = abs(avg_H_zup/avg_H_zincup) - 1
    trans = abs(avg_H_zdown/avg_H_zincdown)
    EB = trans**2 + ref**2     # energy balance
    KD = 1 - EB
    
    return ref, trans, EB, KD