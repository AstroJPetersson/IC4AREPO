import numpy as np

class CollapsingSphere:
	'''
	DOCSTRING
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

	def ambient_medium(self, Nx, Ny, Nz, T, rho):
		# Coordinates:
		x = np.linspace(0, self.boxSize - self.boxSize/Nx, Nx) + self.boxSize/Nx/2
		xx, yy, zz = np.meshgrid(x, x, x)
		Coord = np.zeros((Nx * Ny * Nz, 3))
		Coord[:,0] = xx.flatten()
		Coord[:,1] = yy.flatten()
		Coord[:,2] = zz.flatten()

		# Velocities:
		Vel = np.zeros((Nx * Ny * Nz, 3))

		# Mass:	
		Mass = np.zeros(Nx * Ny * Nz) + (rho * self.boxSize**3 / (Nx * Ny * Nz))
		# Mass = np.full(Nx*Ny*Nz, rho)

		# Internal Energy:
		mu = 1 + 4*0.1
		IntErg = np.full(Nx * Ny * Nz, (3 * self.kb * T) / (2 * mu * self.mp) / self.uinterg)
		
		# IC for ambient medium:
		ic_grid = {
			'Coordinates': Coord,
			'Velocities': Vel,
			'Masses': Mass,
			'InternalEnergy': IntErg,
		}
	
		return ic_grid

	def uniform_sphere(self, N, T, rho, R, J=None):	
		# Coordinates:
		Coord = np.zeros(shape=(N, 3))
		mask = np.full(N, False)
		while len(Coord[mask]) < N:
			Coord[:,0][~mask] = 2*R*np.random.random(N - len(Coord[mask])) - R
			Coord[:,1][~mask] = 2*R*np.random.random(N - len(Coord[mask])) - R
			Coord[:,2][~mask] = 2*R*np.random.random(N - len(Coord[mask])) - R
			mask = (np.linalg.norm(Coord, axis=1) >= 0.1) * (np.linalg.norm(Coord, axis=1) <= R)
		Coord += self.boxSize/2	
	
		# Mass:
		Mass = np.zeros(N) + (rho * (4/3) * np.pi * R**3 / N)
		# Mass = np.full(N, rho)

		# Velocities:
		if J.all() != None:
			# Project positions onto the rotation axis (plane):
			a = Coord - self.boxSize/2
			a1 = (np.dot(a, J) / (np.linalg.norm(J)**2))[:,np.newaxis] * J
			a2 = a - a1
			# Calculate velocities given the total (specific) angular momentum:
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

	def sink(self, x, y, z, vx, vy, vz, m):
		ic_sink = {
			'Coordinates' : np.array([[x, y, z]]),
			'Velocities'  : np.array([[vx, vy, vz]]),
			'Masses'      : np.array([m]),
		}
	
		return ic_sink

	def icgenerate(self, Nx, Ny, Nz, N_sphere, T_grid, T_sphere, rho_grid, rho_sphere, 
				   R_sphere, J_sphere, Coord_sink, Vel_sink, M_sink, prop=False):
		# Generate IC for PartType0:
		ic_grid  = self.ambient_medium(Nx=Nx, Ny=Ny, Nz=Nz, T=T_grid, rho=rho_grid)
		ic_sphere = self.uniform_sphere(N=N_sphere, T=T_sphere, rho=rho_sphere, R=R_sphere, J=J_sphere)

		ic_gas = {}
		for k in ic_grid.keys():
			key_grid = ic_grid[k]
			key_sphere = ic_sphere[k]
			key = np.append(key_grid, key_sphere, axis=0)
			ic_gas[k] = key

		# Generate IC for PartType1:
		ic_sink = self.sink(x=Coord_sink[0], y=Coord_sink[1], z=Coord_sink[2], vx=Vel_sink[0], vy=Vel_sink[1], vz=Vel_sink[2], m=M_sink)

		# Generate ParticleIDs:
		ic_gas['ParticleIDs']  = np.arange(0, N_sphere + Nx * Ny * Nz) + 1
		ic_sink['ParticleIDs'] = np.array([ic_gas['ParticleIDs'][-1] + 1])

		# Collapsing sphere properties:
		if prop == True:
			print('Doing some quick esimates...\n')
			print('For the collapsing sphere:')
			# Sound speed:
			gamma = 5/3
			mu = 1 + 4*0.1
			cs    = ((gamma * self.kb * T_sphere) / (mu * self.mp))**(1/2)

			# Jeans length:
			j = cs / (self.G * rho_sphere * self.udens)**(1/2)

			# Free-fall time:
			tff = ((3 * np.pi) / (32 * self.G * rho_sphere * self.udens))**(1/2)

			print(f'Jeans length [pc] : {j/self.pc}')
			print(f'Sound speed [km s^-1] : {cs/(1e5)}')
			print(f'Free-fall time [Myr]: {tff/(1e6 * 365.25 * 24 * 60 * 60)}')

			print('\nFor the grid cells:')
			# Sound speed:
			gamma = 5/3
			mu = 1 + 4*0.1
			cs    = ((gamma * self.kb * T_grid) / (mu * self.mp))**(1/2)

			# Jeans length:
			j = cs / (self.G * rho_grid * self.udens)**(1/2)

			# Free-fall time:
			tff = ((3 * np.pi) / (32 * self.G * rho_grid * self.udens))**(1/2)

			print(f'Jeans length [pc] : {j/self.pc}')
			print(f'Sound speed [km s^-1] : {cs/(1e5)}')
			print(f'Free-fall time [Myr]: {tff/(1e6 * 365.25 * 24 * 60 * 60)}')

			print('\nChecking generated HDF5-file...\n')
			self.check_ic(f'{self.path}/collapsing_sphere.hdf5')
		
		return ic_gas, ic_sink


