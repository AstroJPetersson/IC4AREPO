import numpy as np
import os
import h5py
import time


def mesh_relaxation(ic, runs, dirjob, icmodify, icmodify_args, filename, savepath, wait='manual'):
    print('  * Now entering the mesh relaxation routine:')
    
    # Move IC-file to mesh relaxation directory:
    os.system(f'mv {ic} {dirjob}')
    
    # Mesh relaxation:
    for i in range(1, runs+1):
        print(f'    - Mesh relaxation: {i}/{runs}')
        
        # -------------- Submission -------------- #
        dirhome = os.getcwd()  # your home directory
        os.chdir(dirjob)       # change directory

        # Submit job:
        os.system('sbatch -Q job.sh')

        # Wait until mesh relaxation is done:
        if wait == 'auto':
            # Automatic waiting:
            print('    - Please Note: You are running mesh relaxation on \'auto\', i.e. the code will automatically\n' +
                  '      check the queue for \'relax\' & remain idle until it no longer runs')
            start = time.time()
            while os.popen('squeue -u jpeterss -h -n relax -o "%.8j"').read():
                time.sleep(1)
                print(f'    - Elapsed time: {int(time.time() - start)} s    ', end='\r')
        else:
            # Manual waiting:
            print('    - Please Note: You are running mesh relaxation on \'manual\',\n' +
                  '      i.e. please keep track of the submitted job manually')
            donejob = 'n'
            while donejob == 'n':
                donejob = input('    - Is the job completed? (y/n): ')
        
        os.chdir(dirhome)  # go back to home directory
        # -------------- Submission -------------- #
        print('    - Mesh relaxation done!')
                
        # Check if all mesh relaxation runs are done:
        if i < runs:
            # If not all runs are done, do some IC-modification(s)
            print('    - Doing some preperations for next mesh relaxation...')
            os.system(f'mv {dirjob}/OUTPUT/snap_001.hdf5 {dirjob}/cloudcollapse.hdf5')
            if icmodify:
                print('    - Doing some modifications of the last snapshot...')
                icmodify(ic=f'{dirjob}/cloudcollapse.hdf5', **icmodify_args)
            print('    - Preperations done!')
        elif i == runs: 
            # If all runs are done, move IC-file to savepath:
            os.system(f'mv {dirjob}/OUTPUT/snap_001.hdf5 {savepath}/{filename}.hdf5')
            with h5py.File(f'{savepath}/{filename}.hdf5', 'r+') as f:
                f['Header'].attrs['Time'] = 0

    print(f'  * All iterations of mesh relaxation done! Congrats!')
    print(f'  * New IC-file: {savepath}/{filename}.hdf5')

    return 0


