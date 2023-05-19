import numpy as np
import os
import h5py
import time


def submit(dirjob, wait='manuel'):
    # Change directory:
    dirhome = os.getcwd()
    os.chdir(dirjob)

    # Submit job:
    os.system('sbatch -Q job.sh')

    if wait == 'auto':
        print('    - Please Note: You are running mesh relaxation on \'auto\', i.e. the code will automatically\n'
               + '      check the queue for \'relax\' & remain idle until it no longer runs')
        while os.popen('squeue -u jpeterss -h -n relax -o "%.8j"').read():
            time.sleep(1)
    else:
        print('    - Please Note: You are running mesh relaxation on \'manual\',\n'
               + '      i.e. please keep track of the submitted job manually')
        donejob = 'n'
        while donejob == 'n':
            donejob = input('    - Is the job completed? (y/n): ')
        
    os.chdir(dirhome)

    return 0

def icmodify_onecloud(ic, centre, R, rho):
    with h5py.File(ic, 'r+') as f:
        # Units:
        ulength = f['Header'].attrs['UnitLength_in_cm']
        umass   = f['Header'].attrs['UnitMass_in_g']
        udens   = umass/(ulength**3)

        # Modify the ICs:
        pos         = f['PartType0']['Coordinates'] * ulength
        mass        = f['PartType0']['Masses'] * umass
        mask        = np.linalg.norm(pos - centre, axis=1) < R
        mass[mask]  = rho
        mass[~mask] = 1e-25
        f['PartType0']['Masses'][:] = mass / udens

    return 0

def mesh_relaxation(ic, runs, dirjob, icmodify, icmodify_args, filename, savepath, wait='manual'):
    print('  * Now entering the mesh relaxation routine:')
    
    os.system(f'mv {ic} {dirjob}')
    
    for i in range(1, runs+1):
        print(f'    - Mesh relaxation: {i}/{runs}')
        submit(dirjob=dirjob, wait=wait)
        print('    - Mesh relaxation done!')

        if i < runs:
            print('    - Doing some preperations for next mesh relaxation...')
            os.system(f'mv {dirjob}/OUTPUT/snap_001.hdf5 {dirjob}/spherecollapse.hdf5')
            if icmodify:
                print('    - Doing some modifications of the last snapshot...')
                icmodify(ic=f'{dirjob}/spherecollapse.hdf5', **icmodify_args)
            print('    - Preperations done!')
        elif i == runs: 
            os.system(f'mv {dirjob}/OUTPUT/snap_001.hdf5 {savepath}/{filename}.hdf5')
            with h5py.File(f'{savepath}/{filename}.hdf5', 'r+') as f:
                f['Header'].attrs['Time'] = 0

    print(f'  * All iterations of mesh relaxation done! Congrats!')
    print(f'  * New IC-file: {savepath}/{filename}.hdf5')

    return 0


