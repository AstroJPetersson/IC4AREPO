import numpy as np
import matplotlib.pyplot as plt
import h5py
from imgcat import imgcat

def check_ic_file(file):
    print('\n\n####### Initializing check_ic_file.py #######\n\n')
    
    # Check generated IC-file:
    with h5py.File(file, 'r') as f:
        print('#### HDF5-file keys: ####')
        for i in f.keys():
            print(i)
    
        print('\n#### HDF5-file header attributes: ####')
        for i in f['Header'].attrs:
            v = f['Header'].attrs[i]
            print(f'{i} : {v}')
    
        print('\n#### PartType0 keys: ####')
        for i in f['PartType0']:
            print(i)
        
        # Data:
        Coord_gas  = f['PartType0']['Coordinates'][:]
        Vel_gas    = f['PartType0']['Velocities'][:]
        Mass_gas   = f['PartType0']['Masses'][:]
        IntErg_gas = f['PartType0']['InternalEnergy'][:]

    # Plot:
    print('\n#### Some simple plots: ####\n')

    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    ax[0].scatter(Coord_gas[:,0], Coord_gas[:,1], marker='.')
    ax[0].set_xlabel('x')
    ax[0].set_ylabel('y')
    ax[0].set_aspect('equal')
    ax[1].scatter(Coord_gas[:,0], Coord_gas[:,2], marker='.')
    ax[1].set_xlabel('x')
    ax[1].set_ylabel('z')
    ax[1].set_aspect('equal')
    ax[2].scatter(Coord_gas[:,1], Coord_gas[:,2], marker='.')
    ax[2].set_xlabel('y')
    ax[2].set_ylabel('z')
    ax[2].set_aspect('equal')
    imgcat(fig)

    fig, ax = plt.subplots(figsize=(6, 6))
    x = np.log10(Mass_gas)
    ax.hist(x, bins=100)
    ax.set_yscale('log')
    ax.set_xlabel('Log10(PartType0 Masses)')
    ax.set_ylabel('# particles')
    imgcat(fig)

    print('\n\n####### Terminating check_ic_file.py #######\n\n')

    return 0


