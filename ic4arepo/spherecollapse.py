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
        Coord = self.boxSize * np.random.rand(N, 3)     

        # Velocities:
        Vel = np.zeros((N, 3))

        # Mass: 
        Mass = np.full(N, rho)

        # Internal Energy:
        mu = 1 + 4*0.1
        IntErg = np.full(N, (3 * self.kb * T) / (2 * mu * self.mp) / self.uinterg)
        
        # IC for ambient medium:
        ic_grid = {
            'Coordinates': Coord,
            'Velocities': Vel,
            'Masses': Mass,
            'InternalEnergy': IntErg,
        }
    
        return ic_grid

    def uniform_sphere(self, N, T, rho, R, J): 
        # Coordinates:
        Coord = R * (2*np.random.rand(N, 3) - 1)        
        mask = np.linalg.norm(Coord, axis=1) > R
        while len(Coord[~mask]) < N:
            Coord[mask] = R * (2*np.random.rand(N - len(Coord[~mask]), 3) - 1)
            mask = np.linalg.norm(Coord, axis=1) > R
        Coord += self.boxSize/2 
    
        # Mass:
        Mass = np.full(N, rho)

        # Velocities:
        if np.shape(J) == 3:
            # Coordinate projection onto the plane of the rotation axis:
            a = Coord - self.boxSize/2
            a1 = (np.dot(a, J) / (np.linalg.norm(J)**2))[:,np.newaxis] * J
            a2 = a - a1
            
            # Velocity of each particle given the total (specific) angular momentum:
            v = (np.linalg.norm(J)/N) * np.linalg.norm(a2, axis=1)
            perp_vec = np.cross(J, a2)
            Vel = v[:,np.newaxis] * perp_vec / (np.linalg.norm(perp_vec, axis=1)[:,np.newaxis])
        
        else:
            Vel = np.zeros((N, 3))
        
        # Internal Energy:
        mu = 1 + 4*0.1
        IntErg = np.full(N, (3 * self.kb * T) / (2 * mu * self.mp) / self.uinterg)
    
        # IC for uniform sphere:
        ic_sphere = {
            'Coordinates': Coord,
            'Velocities': Vel,
            'Masses': Mass,
            'InternalEnergy': IntErg,
        }
        
        return ic_sphere
    
    def sink_particle(self, x, y, z, vx, vy, vz, m):
        ic_sink = {
            'Coordinates' : np.array([[x, y, z]]),
            'Velocities'  : np.array([[vx, vy, vz]]),
            'Masses'      : np.array([m]),
        }
    
        return ic_sink

    def icgenerate(self, N_sphere, T_sphere, rho_sphere, R_sphere, J_sphere, Coord_sink, Vel_sink, M_sink, 
                   N_grid, T_grid, rho_grid, filename, savepath, properties=False, check=False):
        from .write import write_ic_file
        from .check import check_ic_file
        
        # Generate ICs for PartType0:
        ic_sphere = self.uniform_sphere(N=N_sphere, T=T_sphere, rho=rho_sphere, R=R_sphere, J=J_sphere)
        ic_grid  = self.ambient_medium(N=N_grid, T=T_grid, rho=rho_grid)

        ic_gas = {}
        for k in ic_grid.keys():
            key_grid = ic_grid[k]
            key_sphere = ic_sphere[k]
            key = np.append(key_grid, key_sphere, axis=0)
            ic_gas[k] = key

        # Generate ICs for PartType1:
        ic_sink = self.sink_particle(x=Coord_sink[0], y=Coord_sink[1], z=Coord_sink[2], vx=Vel_sink[0], vy=Vel_sink[1], vz=Vel_sink[2], m=M_sink)

        # Generate ParticleIDs:
        ic_gas['ParticleIDs']  = np.arange(0, N_sphere + N_grid) + 1
        ic_sink['ParticleIDs'] = np.array([ic_gas['ParticleIDs'][-1] + 1])
        
        # Write IC-file:
        write_ic_file(filename=filename, savepath=savepath, boxSize=40, partTypes={'PartType0' : ic_gas, 'PartType1' : ic_sink})

        # Check some properties of the collapsing sphere:
        if properties == True:
            print('\n--------------------------------------------------------------------------------')
            print('Quick estimates of the properties of the uniform sphere (and the ambient medium):\n')
            print('Uniform sphere:')
            
            # Sound speed:
            gamma = 5/3
            mu    = 1 + 4*0.1
            cs    = ((gamma * self.kb * T_sphere) / (mu * self.mp))**(1/2)

            # Jeans length:
            j = cs / (self.G * rho_sphere * self.udens)**(1/2)

            # Free-fall time:
            tff = ((3 * np.pi) / (32 * self.G * rho_sphere * self.udens))**(1/2)

            print(f'Jeans length [pc] : {j/self.pc}')
            print(f'Sound speed [km s^-1] : {cs/(1e5)}')
            print(f'Free-fall time [Myr]: {tff/(1e6 * 365.25 * 24 * 60 * 60)}')

            print('\nAmbient medium:')
            
            # Sound speed:
            gamma = 5/3
            mu    = 1 + 4*0.1
            cs    = ((gamma * self.kb * T_grid) / (mu * self.mp))**(1/2)

            # Jeans length:
            j = cs / (self.G * rho_grid * self.udens)**(1/2)

            # Free-fall time:
            tff = ((3 * np.pi) / (32 * self.G * rho_grid * self.udens))**(1/2)

            print(f'Jeans length [pc] : {j/self.pc}')
            print(f'Sound speed [km s^-1] : {cs/(1e5)}')
            print(f'Free-fall time [Myr]: {tff/(1e6 * 365.25 * 24 * 60 * 60)}')
            
            print('--------------------------------------------------------------------------------\n')

        # Check recently generated ic-file:
        if check == True:
            check_ic_file(f'{savepath}/{filename}.hdf5')

        return 0


