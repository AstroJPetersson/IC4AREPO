import numpy as np
import h5py

def write_ic_file(filename, savepath, boxSize, partTypes, massTable=None):
    # Number of particle types:
    nPartTypes = 6

    with h5py.File(f'{savepath + filename}.hdf5', 'w') as f:
        # Write particle types and their respective fields to hdf5-file:
        for ptype in partTypes.keys():
            group = f.create_group(ptype)
            for field in partTypes[ptype]:
                group[field] = partTypes[ptype][field]

        # Set particle counts:
        nPart = np.zeros(nPartTypes, dtype='int64')
        for ptype in partTypes.keys():
            pNum = int(ptype[-1])
            nPart[pNum] = len(partTypes[ptype]['ParticleIDs'])

        # Header:
        h = f.create_group('Header')
        h.attrs['BoxSize'] = boxSize
        h.attrs['NumFilesPerSnapshot'] = 1
        h.attrs['NumPart_ThisFile'] = nPart
        h.attrs['NumPart_Total'] = nPart
        h.attrs['NumPart_Total_HighWord'] = np.zeros(nPartTypes)

        for k in ['Time', 'Redshift', 'Omega0', 'OmegaLambda', 'HubbleParam']:
            h.attrs[k] = 0.0
        for k in ['Sfr', 'Cooling', 'StellarAge', 'Metals', 'Feedback']:
    	    h.attrs[f'Flag_{k}'] = 0
        h.attrs['Flag_DoublePrecision'] = 1

        if massTable is not None:
            h.attrs['MassTable'] = massTable
        else:
            h.attrs['MassTable'] = np.zeros(nPartTypes)

    print(f'  * Initial conditions generated!')
    print(f'  * File name: {filename}')
    print(f'  * Save path: {savepath}')

    return 0


