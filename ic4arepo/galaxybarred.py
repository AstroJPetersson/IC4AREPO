import numpy as np
import matplotlib.pyplot as plt
import h5py
import time
import os
from scipy import interpolate
from imgcat import imgcat

class GalaxyBarred:
    '''
    IC-generator for a barred spiral galaxy
    '''
    def __init__(self, boxSize, ulength, umass, uvel):
        # Simulation set-up:
        self.boxSize  = boxSize
        self.ulength  = ulength
        self.uvel     = uvel
        self.utime    = ulength/uvel
        self.udens    = umass/(ulength**3)
        self.ucoldens = umass/(ulength**2)
        self.uinterg  = uvel**2
        
        # Constants:
        self.kb   = 1.3807e-16     # [erg K^-1]
        self.G    = 6.6726e-8      # [dyne cm^2 g^-2]
        self.mp   = 1.6726e-24     # [g]
        self.pc   = 3.08567758e18  # [cm]
        self.Msol = 1.9891e33      # [g] 
        self.kpc  = 3.08567758e21  # [cm]

        self.iceDir = os.getcwd()

    def galaxy(self, N):
        '''
        Generate artificial ICs 
        '''
        # Coordinates:
        pos = self.boxSize * np.random.rand(N, 3)

        # Velocities: 
        vel = np.zeros((N, 3))
        
        # Densities:
        dens = self.density_profile(pos)
        
        # Internal energy:
        interg = np.zeros(N)

        # Particle IDs:
        iord = np.arange(1, N+1)

        # Initial conditions:
        ic_galaxy = {
            'Coordinates'    : pos,
            'Velocities'     : vel,
            'Masses'         : dens,
            'InternalEnergy' : interg,
            'ParticleIDs'    : iord
        }

        return ic_galaxy

    def velocity_profile(self, pos):
        '''
        Velocity profile for a barred MW-like galaxy
        '''
        R_cut  = 10  * self.kpc                  / self.ulength 
        z_cut  = 150 * self.pc                   / self.ulength
        
        # Load circular velocities from 'check_potential':
        vc = np.loadtxt(f'{self.iceDir}/src/check_potential/vc.txt')

        # Interpolation of ciruclar velocities:
        interp1d = interpolate.interp1d(vc[:,0], vc[:,1])
        vel_c = interp1d(np.linalg.norm(pos, axis=1))
        
        # Cut:
        x, y, z = pos.T - self.boxSize/2
        R = np.sqrt(x**2 + y**2)
        vel_c[(R > R_cut) * (np.abs(z) > z_cut)] = 0
        
        # Velocity vector:
        jaxis = np.array([0, 0, -1])
        pos_proj = (pos - self.boxSize/2) - np.dot(pos, jaxis)[:,np.newaxis] * jaxis
        vel_t = np.cross(jaxis, pos_proj)
        vel = vel_c[:,np.newaxis] * (vel_t / (np.linalg.norm(vel_t, axis=1)[:,np.newaxis]))

        return vel

    def density_profile(self, pos):
        '''
        Density profile for a barred MW-like galaxy (ref: Robin Tress) 
        '''
        Sigma0 = 50  * self.Msol/self.pc/self.pc / self.ucoldens
        zd     = 85  * self.pc                   / self.ulength
        Rm     = 1.5 * self.kpc                  / self.ulength
        Rd     = 7   * self.kpc                  / self.ulength 
        R_cut  = 10  * self.kpc                  / self.ulength
        z_cut  = 150 * self.pc                   / self.ulength

        x, y, z = pos.T - self.boxSize/2
        R = np.sqrt(x**2 + y**2)
        density = Sigma0/(4*zd) * np.exp(-Rm/R - R/Rd) / np.cosh(z/(2*zd))**2
        density[R > R_cut] = 1e-28 / self.udens
        density[np.abs(z) > z_cut] = 1e-28 / self.udens

        return density
    
    def meshrelax(self, jobDir, wait='manuel', sleep=120):
        '''
        Mesh-relaxation of IC-file
        '''
        os.chdir(jobDir)

        print('Submitting and waiting for job to be completed...')
        os.system('sbatch job.sh')
        
        # 'sleep' wait method:
        if wait == 'sleep':
            if os.path.isfile('Arepo.out'):
                with open('Arepo.out') as A:
                    a = A.readlines()
                    line_cut = a[-4][13:-10]
                    last_run = float(line_cut)
                    sleep = 1.3*last_run
                    print(f'Last code run: {last_run} s | Estimated code run: {sleep} s')
            time.sleep(sleep)
        # 'auto' wait method:
        elif wait == 'auto':
            print('You are running with wait=\'auto\', the code will now automatically check the ' +
                  'slurm queue for the job \'relax\' and remain idle until the job no longer runs')
            while os.popen('squeue -u jpeterss -h -n relax -o "%.8j"').read():
                time.sleep(1)
        # 'manual' wait method:
        else:
            print('(please keep track of the submitted job manually)')
            proc = 'n'
            while proc == 'n':
                proc = input('Is the job completed? (y/n): ')
        
        os.chdir(self.iceDir)

        return 0

    def icmodify(self, jobDir):
        '''
        Modification of the final snapshot in mesh-relaxation
        '''
        # Move final snapshot from previous mesh-relax and rename it to galaxy.hdf5
        os.system(f'mv {jobDir}/OUTPUT/snap_001.hdf5 {jobDir}/galaxy_barred.hdf5')

        # Modify the new IC-file with the density profile:
        with h5py.File(f'{jobDir}/galaxy_barred.hdf5', 'r+') as f:
            pos = f['PartType0']['Coordinates'][:]
            dens_upd = self.density_profile(pos)
            f['PartType0']['Masses'][:] = dens_upd

        return 0

    def icgenerate(self, N, runs, jobDir, saveDir, wait='manual', sleep=120):
        '''
        Generate IC-file 
        '''     
        from .write import write_ic_file

        # Generate galaxy model and wrtie ic-file:
        print('Generating artifical ICs for a barred spiral Milky Way-like galaxy...')
        ic_galaxy = self.galaxy(N)
        write_ic_file(filename='galaxy_barred', savepath=jobDir, boxSize=self.boxSize, partTypes={'PartType0' : ic_galaxy})

        # Run and repeat mesh-relaxation:
        print(f'\nInitialising mesh-relaxtion on the generated IC-file:')
        for i in range(1, runs+1):
            print(f'Mesh-relaxation {i}/{runs}: Running...')
            self.meshrelax(jobDir=jobDir, wait=wait, sleep=sleep)
            print(f'Mesh-relaxation {i}/{runs}: DONE!\n')

            if i < runs:
                print('Preparing IC-file for next mesh-relaxation...')
                self.icmodify(jobDir=jobDir)
                print('Preparation done!\n')
            elif i == runs: 
                os.system(f'mv {jobDir}/OUTPUT/snap_001.hdf5 {saveDir}/galaxy_barred.hdf5')
                
                # Update the final IC-file with the ciruclar velocity of each particle:
                with h5py.File(f'{saveDir}/galaxy_barred.hdf5', 'r+') as f:
                    pos = f['PartType0']['Coordinates'][:]
                    vel_upd = self.velocity_profile(pos)
                    f['PartType0']['Velocities'][:] = vel_upd

                print(f'All interations of mesh-relaxation done! New IC-file \'galaxy_barred\' can be found at: \'{saveDir}\'')

        return 0


