import numpy as np
import os
import h5py

from .write import write_ic_file
from .check import check_ic_file
from .relax import mesh_relaxation


class Streams:
    '''
    Streams: initial conditions for two opposite streams
    '''
    def __init__(self, boxSize, ulength, umass, uvel):
        # Simulation set-up:
        self.boxSize  = boxSize
        self.ulength  = ulength
        self.umass    = umass
        self.uvel     = uvel
        self.utime    = ulength/uvel
        self.udens    = umass/(ulength**3)
        self.ucoldens = umass/(ulength**2)
        self.uinterg  = uvel**2
        self.uangmom  = ulength * uvel

        # Constants:
        self.kb   = 1.3807e-16     # [erg K^-1]
        self.G    = 6.6726e-8      # [dyne cm^2 g^-2]
        self.mp   = 1.6726e-24     # [g]
        self.pc   = 3.08567758e18  # [cm]
        self.Msol = 1.9891e33      # [g]
        self.kpc  = 3.08567758e21  # [cm]
    
    
    #-------------------------------------------------------------------------#
    def streams(self, N, rho, T, v):
         # Coordinates:
        pos = self.boxSize * np.random.rand(N, 3)

        # Masses:
        mass = np.full(N, rho)

        # Velocities:
        vel = self.velocity(pos, v)

        # Internal Energy:
        mu = 1 + 4*0.1
        interg = np.full(N, (3 * self.kb * T) / (2 * mu * self.mp))

        ic_stream = {
            'Coordinates': pos / self.ulength,
            'Velocities': vel / self.uvel,
            'Masses': mass / self.udens,
            'InternalEnergy': interg / self.uinterg,
        }
        
        return ic_streams
    
    
    def velocity(self, pos, v):
        vel = np.zeros((N, 3))
        stream = np.array([v, 0, 0])
        
        mask       = pos[:,1] > self.boxSize / 2
        vel[mask]  = -stream
        vel[~mask] = stream
        
        return vel
    
    
    def icmodify(rho):
        # Modify the ICs:
        with h5py.File(ic, 'r+') as f:
            f['PartType0']['Masses'][:] = np.full(len(mass), rho) / self.udens
        
        return None
    
    
    def sink_particle(self, x, y, z, vx, vy, vz, m):
        ic_sink = {
            'Coordinates' : np.array([[x, y, z]]) / self.ulength,
            'Velocities'  : np.array([[vx, vy, vz]]) / self.uvel,
            'Masses'      : np.array([m]) / self.umass,
         }
    
        return ic_sink
    
    
    #-------------------------------------------------------------------------#
    def generate(self, N, rho, T, vel_gas, pos_sink, vel_sink, mass_sink, filename, savepath, check=False, relax=False, N_relax=1, wait='manual'):
        # Generate ICs for PartType0:
        ic_gas = self.streams(N=N, rho=rho, T=T, v=vel_gas)

        # Generate ICs for PartType5:
        ic_sink = self.sink_particle(x=pos_sink[0], y=pos_sink[1], z=pos_sink[2], vx=vel_sink[0], vy=vel_sink[1], 
                                     vz=vel_sink[2], m=mass_sink)

        # Generate ParticleIDs:
        ic_gas['ParticleIDs']  = np.arange(0, N) + 1
        ic_sink['ParticleIDs'] = np.array([ic_gas['ParticleIDs'][-1] + 1])
        
        # Write IC-file:
        write_ic_file(filename=filename, savepath=savepath, boxSize=self.boxSize / self.ulength, 
                      partTypes={'PartType0' : ic_gas, 'PartType5' : ic_sink})

        # Check IC-file:
        if check == True:
            check_ic_file(f'{savepath}/{filename}.hdf5')

        # Mesh relaxation:
        if relax == True:
            dirjob        = os.getcwd() + '/src/mesh_relaxation/cloudcollapse/'
            icmodify      = self.icmodify()
            icmodify_args = {'rho' : rho}
            mesh_relaxation(ic=f'{savepath}/{filename}.hdf5', runs=N_relax, dirjob=dirjob, icmodify=icmodify, 
                            icmodify_args=icmodify_args, filename=filename, savepath=savepath, wait=wait)
           
            # Update of the velocity:
            with h5py.File(f'{savepath}/{filename}.hdf5', 'r+') as f:
                pos = f['PartType0']['Coordinates'] * ulength
                vel = self.velocity(pos=pos, v=vel_gas)
                f['PartType0']['Velocities'][:] = vel / uvel
        
        return None
    
    