import numpy as np
import matplotlib.pyplot as plt
import h5py
from imgcat import imgcat

def check_ic_file(file):
    # Check generated IC-file:
    with h5py.File(file, 'r') as f:
        print('HDF5-file keys:')
        for i in f.keys():
            print(i)
    
        print('\nHDF5-file header attributes:')
        for i in f['Header'].attrs:
            v = f['Header'].attrs[i]
            print(f'{i} : {v}')
    
        print('\nPartType0 keys:')
        for i in f['PartType0']:
            print(i)
        
        # Data:
        Coord_gas  = f['PartType0']['Coordinates'][:]
        Vel_gas    = f['PartType0']['Velocities'][:]
        Mass_gas   = f['PartType0']['Masses'][:]
        IntErg_gas = f['PartType0']['InternalEnergy'][:]

    print(f'Max Mass: {np.max(Mass_gas)}')
    print(f'Min Mass: {np.min(Mass_gas)}')

    # Plot:
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.hist(np.log10(Mass_gas), bins=100)
    ax.set_yscale('log')
    imgcat(fig)

    '''
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.hist(np.linalg.norm(Coord_gas[:,0:3], axis=1), bins=100)
    ax.set_yscale('log')
    ax.set_xlim(0, 500)
    imgcat(fig)
    '''

    return 0


