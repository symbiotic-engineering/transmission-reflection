def wave_height(total, incoming_fse, lam, xtrans, ytrans, farm, rel_dim, nx, ny, x1, x2, y1, y2):
    import numpy as np
    import matplotlib.pyplot as plt
    import warnings
    warnings.filterwarnings("ignore", category=np.ComplexWarning)
    
    ##################################################################
    # Extract relevant columns
    mid_y = int(ny / 2)
    mid_x = int(nx / 2)
    convx = int(nx/(abs(x1)+x2))      # to convert meters to grid points
    convy = int(ny/(abs(y1)+y2))
    
    # Extract z_up and z_down using rel_dim to determine the region
    zinc_up = incoming_fse[mid_y, :mid_x - int(rel_dim*convx)]
    zinc_down = incoming_fse[mid_y, mid_x + int(rel_dim*convx):]
    z_up = total[mid_y, :mid_x - int(rel_dim*convx)]
    z_down = total[mid_y, mid_x + int(rel_dim*convx):]
    
    avg_H_zup = np.array([np.mean(abs(z_up))])
    avg_H_zincup = np.array([np.mean(abs(zinc_up))])
    avg_H_zdown = np.array([np.mean(abs(z_down))])
    avg_H_zincdown = np.array([np.mean(abs(zinc_down))])
    
    if farm:
        z_upWEC1 = total[mid_y + int(ytrans[0]*convy), :mid_x - int(rel_dim*convx) + int(xtrans[0]*convx)]
        z_upWEC3 = total[mid_y + int(ytrans[1]*convy), :mid_x - int(rel_dim*convx) + int(xtrans[1]*convx)]

        z_downWEC1 = total[mid_y + int(ytrans[0]*convy), mid_x + int(rel_dim*convx) + int(xtrans[0]*convx):]
        z_downWEC3 = total[mid_y + int(ytrans[1]*convy), mid_x + int(rel_dim*convx) + int(xtrans[1]*convx):]

        avg_H_zup = np.array([np.mean(abs(z_upWEC1)), np.mean(abs(z_up)), np.mean(abs(z_upWEC3))])
        avg_H_zdown = np.array([np.mean(abs(z_downWEC1)), np.mean(abs(z_down)), np.mean(abs(z_downWEC3))])

    ref = abs(avg_H_zup / avg_H_zincup) - 1
    trans = abs(avg_H_zdown / avg_H_zincdown)
    EB = trans**2 + ref**2     # energy balance
    KD = 1 - EB

    # ##################################################################
    # # Plotting z_up and z_down geographically

    # x = np.linspace(x1, x2, nx)
    # y = np.linspace(y1, y2, ny)
    # X, Y = np.meshgrid(x, y)
    
    # plt.figure(figsize=(12, 8))
    
    # # Contour plot of incoming_fse
    # contour = plt.contourf(X, Y, incoming_fse, cmap='viridis')
    # plt.colorbar(contour, label='Incoming FSE')
    
    # # Plot z_up line
    # x_up = np.linspace(x1, x1 + (mid_x - (rel_dim*convx)) * (x2 - x1) / nx, mid_x - (rel_dim*convx))
    # plt.plot(x_up, [0] * len(x_up), label='z_up', color='blue')
    
    # # Plot z_down line
    # x_down = np.linspace(x1 + (mid_x + (rel_dim*convx)) * (x2 - x1) / nx, x2, nx - mid_x - (rel_dim*convx))
    # plt.plot(x_down, [0] * len(x_down), label='z_down', color='red')
    
    # if farm:
    #     # Calculate new positions for WEC lines based on xtrans and ytrans
    #     y_upWEC1 = mid_y + int(ytrans[0]*convy)
    #     y_upWEC3 = mid_y + int(ytrans[1]*convy)
        
    #     x_upWEC1 = np.linspace(x1, x1 + (mid_x - (rel_dim*convx) + int(xtrans[0]*convx)) * (x2 - x1) / nx, mid_x - (rel_dim*convx) + int(xtrans[0]*convx))
    #     x_upWEC3 = np.linspace(x1, x1 + (mid_x - (rel_dim*convx) + int(xtrans[1]*convx)) * (x2 - x1) / nx, mid_x - (rel_dim*convx) + int(xtrans[1]*convx))

    #     y_downWEC1 = mid_y + int(ytrans[0]*convy)
    #     y_downWEC3 = mid_y + int(ytrans[1]*convy)
        
    #     x_downWEC1 = np.linspace(x1 + (mid_x + (rel_dim*convx) + int(xtrans[0]*convx)) * (x2 - x1) / nx, x2, nx - mid_x - (rel_dim*convx) - int(xtrans[0]*convy))
    #     x_downWEC3 = np.linspace(x1 + (mid_x + (rel_dim*convx) + int(xtrans[1]*convx)) * (x2 - x1) / nx, x2, nx - mid_x - (rel_dim*convx) - int(xtrans[1]*convy))
        
    #     plt.plot(x_upWEC1, [y_upWEC1 * (y2 - y1) / ny + y1] * len(x_upWEC1), 'b--', label='z_upWEC1')
    #     plt.plot(x_upWEC3, [y_upWEC3 * (y2 - y1) / ny + y1] * len(x_upWEC3), 'b-.', label='z_upWEC3')
    #     plt.plot(x_downWEC1, [y_downWEC1 * (y2 - y1) / ny + y1] * len(x_downWEC1), 'r--', label='z_downWEC1')
    #     plt.plot(x_downWEC3, [y_downWEC3 * (y2 - y1) / ny + y1] * len(x_downWEC3), 'r-.', label='z_downWEC3')
    
    # plt.xlabel('X Coordinate')
    # plt.ylabel('Y Coordinate')
    # plt.title('Geographical Plot of z_up and z_down on Incoming FSE Contour')
    # plt.legend()
    # plt.grid(True)
    # plt.show()
    # print('wavesss')
    # plt.savefig('locations.pdf')
    # print('woesss')
    
    return ref, trans, EB, KD

# def wave_height(total, incoming_fse,lam,xtrans,ytrans,farm,rel_dim,nx,ny,x1,x2,y1,y2):
#     import numpy as np
#     from scipy.signal import find_peaks
#     import pandas as pd
#     import matplotlib.pyplot as plt
#     import warnings
#     warnings.filterwarnings("ignore", category=np.ComplexWarning)
#     ##################################################################
#     # Extract relevant columns
#     zinc_up = incoming_fse[int(ny/2),:int((nx/2)-rel_dim)]
#     zinc_down = incoming_fse[int(ny/2),int((nx/2)+rel_dim):]
#     z_up = total[int(ny/2),:int((nx/2)-rel_dim)]
#     z_down = total[int(ny/2),int((nx/2)+rel_dim):]

#     avg_H_zup = np.array([np.mean(abs(z_up))])
#     avg_H_zincup = np.array([np.mean(abs(zinc_up))])
#     avg_H_zdown = np.array([np.mean(abs(z_down))])
#     avg_H_zincdown = np.array([np.mean(abs(zinc_down))])


#     if farm == True:
#         z_upWEC1 = total[int((ny/2) + (ytrans[0]/2)),:int((nx/2)-rel_dim+(xtrans[0]/2))]
#         z_upWEC3 = total[int((ny/2) + (ytrans[1]/2)),:int((nx/2)-rel_dim+(xtrans[1]/2))]

#         z_downWEC1 = total[int((ny/2) + (ytrans[0]/2)),int((nx/2)+rel_dim+(xtrans[0]/2)):]
#         z_downWEC3 = total[int((ny/2) + (ytrans[1]/2)),int((nx/2)+rel_dim+(xtrans[1]/2)):]

#         avg_H_zup = np.array([np.mean(abs(z_upWEC1)),np.mean(abs(z_up)),np.mean(abs(z_upWEC3))])
#         avg_H_zdown = np.array([np.mean(abs(z_downWEC1)),np.mean(abs(z_down)),np.mean(abs(z_downWEC3))])

#     ref = abs(avg_H_zup/avg_H_zincup) - 1
#     trans = abs(avg_H_zdown/avg_H_zincdown)
#     EB = trans**2 + ref**2     # energy balance
#     KD = 1 - EB

#     # Creating coordinates for z_up and z_down
#     z_up_x = np.linspace(x1, x2, int((nx/2)-rel_dim))
#     z_up_y = np.full_like(z_up_x, y1)
#     z_down_x = np.linspace(x1, x2, int((nx/2)-rel_dim))
#     z_down_y = np.full_like(z_down_x, y2)
#     if farm:
#         # Add coordinates for farm locations
#         z_up_x = np.concatenate([z_up_x, xtrans])
#         z_up_y = np.concatenate([z_up_y, ytrans])
#         z_down_x = np.concatenate([z_down_x, xtrans])
#         z_down_y = np.concatenate([z_down_y, ytrans])
    
#     plt.figure(figsize=(10, 6))
    
#     # Plot z_up points
#     plt.scatter(z_up_x, z_up_y, c='blue', label='z_up', alpha=0.6)
    
#     # Plot z_down points
#     plt.scatter(z_down_x, z_down_y, c='red', label='z_down', alpha=0.6)
    
#     plt.xlabel('X Coordinate')
#     plt.ylabel('Y Coordinate')
#     plt.title('Geographical Plot of z_up and z_down')
#     plt.legend()
#     plt.grid(True)
#     print('wavesss')
#     plt.savefig('locations.pdf')
#     print('woesss')
    
#     return ref, trans, EB, KD