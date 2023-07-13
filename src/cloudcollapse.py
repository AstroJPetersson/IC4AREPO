import numpy as np
import os
import h5py

from .write import write_ic_file
from .check import check_ic_file
from .relax import mesh_relaxation


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


def icmodify_r2cloud(ic, centre, R, rho):
    pc = 3.08567758e18  #[cm]
    
    with h5py.File(ic, 'r+') as f:
        # Units:
        ulength = f['Header'].attrs['UnitLength_in_cm']
        umass   = f['Header'].attrs['UnitMass_in_g']
        udens   = umass/(ulength**3)

        # Modify the ICs:
        pos          = f['PartType0']['Coordinates'] * ulength
        mass         = f['PartType0']['Masses'] * umass
        radius       = np.linalg.norm(pos - centre, axis=1)
        mask1        = radius < R
        mask2        = radius < 1.0 * pc
        mass[mask1]  = rho * (1.0 * pc / radius[mask1])**2 
        mass[mask2]  = rho
        mass[~mask1] = 1e-28
        f['PartType0']['Masses'][:] = mass / udens

    return 0


class CloudCollapse:
    '''
    CloudCollapse: initial conditions for collapsing cloud(s)
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
    def uniform_cloud(self, N, centre, kick, R, rho, T, P, rot):
        # Coordinates:
        pos = self.boxSize * np.random.rand(N, 3)

        # Mass:
        mass = np.full(N, 1e-25)
        mask = np.linalg.norm(pos - centre, axis=1) < R
        mass[mask] = rho

        # Velocities:
        vel = self.velocity(pos - centre, R, kick, P, rot)

        # Internal Energy:
        mu = 1 + 4*0.1
        interg = np.full(N, (3 * self.kb * 1e4) / (2 * mu * self.mp))
        interg[mask] = (3 * self.kb * T) / (2 * mu * self.mp)

        ic_cloud = {
            'Coordinates': pos / self.ulength,
            'Velocities': vel / self.uvel,
            'Masses': mass / self.udens,
            'InternalEnergy': interg / self.uinterg,
        }

        return ic_cloud

    def velocity(self, pos, R, kick, P, rot):
        # Velocity:
        if P > 0:
            # Rotation axis:
            rot = rot / np.linalg.norm(rot)

            # Coordinate projection onto the plane of the rotation axis:
            pos1 = np.dot(pos, rot)[:,np.newaxis] * rot
            pos2 = pos - pos1

            # Velocity vectors based on the angular velocity:
            vec  = np.cross(rot, pos2)
            vec  = vec / np.linalg.norm(vec)
            w    = 2 * np.pi / P
            v    = w * np.linalg.norm(pos2, axis=1)
            vel  = v[:,np.newaxis] * vec
            mask = np.linalg.norm(pos, axis=1) < R
            vel[~mask] = 0
        else:
            vel = np.zeros(np.shape(pos))
        
        vel += kick

        return vel
    
    #-------------------------------------------------------------------------#
    def r2_cloud(self, N, centre, kick, R, rho, T):
        # Coordinates:
        N_cloud    = int(0.9 * N)
        N_medium   = int(0.1 * N)
        pos_cloud  = R * (2 * np.random.rand(N_cloud, 3) - 1) + centre
        pos_medium = self.boxSize * np.random.rand(N_medium, 3)
        pos        = np.append(pos_cloud, pos_medium, axis=0)

        # Mass:
        dens        = np.full(N, 1e-28)
        radius      = np.linalg.norm(pos - centre, axis=1)
        mask1       = radius < R
        mask2       = radius < 1.0 * self.pc
        dens[mask1] = rho * (1.0 * self.pc / radius[mask1])**2 
        dens[mask2] = rho

        # Velocities:
        vel  = np.zeros(np.shape(pos))
        vel += kick

        # Internal Energy:
        mu = 1 + 4*0.1
        interg = np.full(N, (3 * self.kb * 1e4) / (2 * mu * self.mp))
        interg[mask1] = (3 * self.kb * T) / (2 * mu * self.mp)

        ic_cloud = {
            'Coordinates': pos / self.ulength,
            'Velocities': vel / self.uvel,
            'Masses': dens / self.udens,
            'InternalEnergy': interg / self.uinterg,
        }

        return ic_cloud

    #-------------------------------------------------------------------------#
    def two_uniform_clouds(self, N, centre1, centre2, kick1, kick2, R1, R2, rho1, rho2, T1, T2):
        # Coordinates:
        pos = self.boxSize * np.random.rand(N, 3)

        # Mass:
        mass = np.full(N, 1e-24)
        mask1 = np.linalg.norm(pos - centre1, axis=1) < R1
        mask2 = np.linalg.norm(pos - centre2, axis=1) < R2
        mass[mask1] = rho1
        mass[mask2] = rho2
        
        # Velocities:
        vel = np.zeros((N, 3))
        vel[mask1] += kick1
        vel[mask2] += kick2

        # Internal Energy:
        mu = 1 + 4*0.1
        interg = np.full(N, (3 * self.kb * 1e4) / (2 * mu * self.mp))
        interg[mask1] = (3 * self.kb * T1) / (2 * mu * self.mp)
        interg[mask2] = (3 * self.kb * T2) / (2 * mu * self.mp)

        ic_clouds = {
            'Coordinates': pos / self.ulength,
            'Velocities': vel / self.uvel,
            'Masses': mass / self.udens,
            'InternalEnergy': interg / self.uinterg,
        }
        
        return ic_clouds

    #-------------------------------------------------------------------------#
    def sink_particle(self, x, y, z, vx, vy, vz, m):
        ic_sink = {
            'Coordinates' : np.array([[x, y, z]]) / self.ulength,
            'Velocities'  : np.array([[vx, vy, vz]]) / self.uvel,
            'Masses'      : np.array([m]) / self.umass,
         }
    
        return ic_sink

    #-------------------------------------------------------------------------#
    def icgenerate_uniform_cloud(self, N, pos_cloud, vel_cloud, R, rho, T, P, rot, pos_sink, vel_sink, mass_sink, 
                                 filename, savepath, check=False, relax=False, N_relax=1, wait='manual'):
        # Generate ICs for PartType0:
        ic_gas = self.uniform_cloud(N=N, centre=pos_cloud, kick=vel_cloud, R=R, rho=rho, T=T, P=P, rot=rot)

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
            icmodify      = icmodify_onecloud
            icmodify_args = {'centre' : pos_cloud, 'R' : R, 'rho' : rho}
            mesh_relaxation(ic=f'{savepath}/{filename}.hdf5', runs=N_relax, dirjob=dirjob, icmodify=icmodify, 
                            icmodify_args=icmodify_args, filename=filename, savepath=savepath, wait=wait)
           
            # Final update of the velocity:
            with h5py.File(f'{savepath}/{filename}.hdf5', 'r+') as f:
                ulength = f['Header'].attrs['UnitLength_in_cm']
                uvel    = f['Header'].attrs['UnitVelocity_in_cm_per_s']
                pos     = f['PartType0']['Coordinates'] * ulength 
                pos     = pos - pos_cloud 
                vel     = self.velocity(pos=pos, R=R, kick=vel_cloud, P=P, rot=rot)
                f['PartType0']['Velocities'][:] = vel / uvel

        return None
    
    #-------------------------------------------------------------------------#
    def icgenerate_r2_cloud(self, N, pos_cloud, vel_cloud, R, rho, T, pos_sink, vel_sink, mass_sink, 
                            filename, savepath, check=False, relax=False, N_relax=1, wait='manual'):
        # Generate ICs for PartType0:
        ic_gas = self.r2_cloud(N=N, centre=pos_cloud, kick=vel_cloud, R=R, rho=rho, T=T)

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
            icmodify      = icmodify_r2cloud
            icmodify_args = {'centre' : pos_cloud, 'R' : R, 'rho' : rho}
            mesh_relaxation(ic=f'{savepath}/{filename}.hdf5', runs=N_relax, dirjob=dirjob, icmodify=icmodify, 
                            icmodify_args=icmodify_args, filename=filename, savepath=savepath, wait=wait)
           
            # Final update of the velocity:
            with h5py.File(f'{savepath}/{filename}.hdf5', 'r+') as f:
                vel = f['PartType0']['Velocities'][:]
                vel = np.zeros(np.shape(vel))
                f['PartType0']['Velocities'][:] = vel
        
        return None

    #-------------------------------------------------------------------------#
    def icgenerate_two_uniform_clouds(self, N, pos_cloud1, pos_cloud2, vel_cloud1, vel_cloud2, R1, R2, rho1, rho2, T1, T2, 
                                      pos_sink, vel_sink, mass_sink, filename, savepath, check=False, relax=False, 
                                      N_relax=1, wait='manual'):
        # Generate ICs for PartType0:
        ic_gas = self.twoclouds(N=N, centre1=pos_cloud1, centre2=pos_cloud2, kick1=vel_cloud1, kick2=vel_cloud2, 
                                R1=R1, R2=R2, rho1=rho1, rho2=rho2, T1=T1, T2=T2)

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
            icmodify      = None
            icmodify_args = None
            mesh_relaxation(ic=f'{savepath}/{filename}.hdf5', runs=N_relax, dirjob=dirjob, icmodify=icmodify, 
                            icmodify_args=icmodify_args, filename=filename, savepath=savepath, wait=wait)

        return None


