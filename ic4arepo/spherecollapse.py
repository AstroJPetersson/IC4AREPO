import numpy as np

class SphereCollapse:
    '''
    SphereCollapse
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

    def ambient_medium(self, N, T, rho):
        # Coordinates:
        pos = self.boxSize * np.random.rand(N, 3)     

        # Velocities:
        vel = np.zeros((N, 3))

        # Mass: 
        mass = np.full(N, rho)

        # Internal Energy:
        mu = 1 + 4*0.1
        interg = np.full(N, (3 * self.kb * T) / (2 * mu * self.mp) / self.uinterg)
        
        # IC for ambient medium:
        ic_grid = {
            'Coordinates': pos / self.ulength,
            'Velocities': vel / self.uvel,
            'Masses': mass / self.udens,
            'InternalEnergy': interg / self.uinterg,
        }
    
        return ic_grid

    def uniform_sphere(self, N, T, rho, R, J): 
        # Coordinates:
        pos = R * (2*np.random.rand(N, 3) - 1)        
        mask = np.linalg.norm(pos, axis=1) > R
        while len(pos[~mask]) < N:
            pos[mask] = R * (2*np.random.rand(N - len(pos[~mask]), 3) - 1)
            mask = np.linalg.norm(pos, axis=1) > R
        pos += self.boxSize/2 
    
        # Mass:
        mass = np.full(N, rho)

        # Velocities:
        if np.shape(J) == 3:
            # Coordinate projection onto the plane of the rotation axis:
            a = pos - self.boxSize/2
            a1 = (np.dot(a, J) / (np.linalg.norm(J)**2))[:,np.newaxis] * J
            a2 = a - a1
            
            # Velocity of each particle given the total (specific) angular momentum:
            v = (np.linalg.norm(J)/N) * np.linalg.norm(a2, axis=1)
            perp_vec = np.cross(J, a2)
            vel = v[:,np.newaxis] * perp_vec / (np.linalg.norm(perp_vec, axis=1)[:,np.newaxis])
        else:
            vel = np.zeros((N, 3))
        
        # Internal Energy:
        mu = 1 + 4*0.1
        interg = np.full(N, (3 * self.kb * T) / (2 * mu * self.mp) / self.uinterg)
    
        # IC for uniform sphere:
        ic_sphere = {
            'Coordinates': pos / self.ulength,
            'Velocities': vel / self.uvel,
            'Masses': mass / self.udens,
            'InternalEnergy': interg / self.uinterg,
        }
        
        return ic_sphere
    
    def sink_particle(self, x, y, z, vx, vy, vz, m):
        ic_sink = {
            'Coordinates' : np.array([[x, y, z]]) / self.ulength,
            'Velocities'  : np.array([[vx, vy, vz]]) / self.uvel,
            'Masses'      : np.array([m]) / self.umass,
         }
    
        return ic_sink

    def icgenerate(self, N_sphere, T_sphere, rho_sphere, R_sphere, J_sphere, pos_sink, vel_sink, mass_sink, 
                   N_grid, T_grid, rho_grid, filename, savepath, jeans_check=False, check=False):
        from .write import write_ic_file
        from .check import check_ic_file
        
        # Generate ICs for PartType0:
        ic_sphere = self.uniform_sphere(N=N_sphere, T=T_sphere, rho=rho_sphere, R=R_sphere, J=J_sphere)
        ic_grid   = self.ambient_medium(N=N_grid, T=T_grid, rho=rho_grid)

        ic_gas = {}
        for k in ic_grid.keys():
            key_grid = ic_grid[k]
            key_sphere = ic_sphere[k]
            key = np.append(key_grid, key_sphere, axis=0)
            ic_gas[k] = key

        # Generate ICs for PartType5:
        ic_sink = self.sink_particle(x=pos_sink[0], y=pos_sink[1], z=pos_sink[2], vx=vel_sink[0], vy=vel_sink[1], vz=vel_sink[2], m=mass_sink)

        # Generate ParticleIDs:
        ic_gas['ParticleIDs']  = np.arange(0, N_sphere + N_grid) + 1
        ic_sink['ParticleIDs'] = np.array([ic_gas['ParticleIDs'][-1] + 1])
        
        # Write IC-file:
        write_ic_file(filename=filename, savepath=savepath, boxSize=40, partTypes={'PartType0' : ic_gas, 'PartType5' : ic_sink})

        # Do some checks on the Jeans length of the collapsing sphere & ambient medium:
        if jeans_check == True:
            print(f'  * Since jeans_check=True, do some quick estimates of Jeans length of the uniform sphere (& ambient medium)')
            print('  * Uniform sphere:')
            
            # Sound speed:
            gamma = 5/3
            mu    = 1 + 4*0.1
            cs    = ((gamma * self.kb * T_sphere) / (mu * self.mp))**(1/2)

            # Jeans length:
            j = cs / (self.G * rho_sphere)**(1/2)

            # Free-fall time:
            tff = ((3 * np.pi) / (32 * self.G * rho_sphere))**(1/2)

            print(f'    Jeans length [pc] : {j/self.pc}')
            print(f'    Sound speed [km s^-1] : {cs/(1e5)}')
            print(f'    Free-fall time [Myr]: {tff/(1e6 * 365.25 * 24 * 60 * 60)}')

            print('  * Ambient medium:')
            
            # Sound speed:
            gamma = 5/3
            mu    = 1 + 4*0.1
            cs    = ((gamma * self.kb * T_grid) / (mu * self.mp))**(1/2)

            # Jeans length:
            j = cs / (self.G * rho_grid)**(1/2)

            # Free-fall time:
            tff = ((3 * np.pi) / (32 * self.G * rho_grid))**(1/2)

            print(f'    Jeans length [pc] : {j/self.pc}')
            print(f'    Sound speed [km s^-1] : {cs/(1e5)}')
            print(f'    Free-fall time [Myr]: {tff/(1e6 * 365.25 * 24 * 60 * 60)}')

        # Check IC-file:
        if check == True:
            check_ic_file(f'{savepath}/{filename}.hdf5')

        return 0


