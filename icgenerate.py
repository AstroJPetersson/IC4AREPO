#!/home/astro/jpeterss/anaconda3/bin/python

# -------------- Packages
import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from imgcat import imgcat

# ------------- IC generate class
class IceTap:
	'''
	ICETAP : Initial Conditions gEneraTor for Arepo
	'''
	def __init__(self, path, boxsize, ulength, umass, uvel):
		# Simulation set-up:
		self.path    = path
		self.boxsize = boxsize
		self.ulength = ulength
		self.umass   = umass
		self.uvel    = uvel
		self.utime   = ulength/uvel
		self.udens   = umass/(ulength**3)
		self.uinterg = uvel**2
		self.uangmom = ulength * uvel

		# Constants:
		self.kb = 1.3807e-16     # [erg K^-1]
		self.G  = 6.6726e-8      # [dyne cm^2 g^-2]
		self.mp = 1.6726e-24     # [g]
		self.pc = 3.08567758e18  # [cm]

	def ambient_medium(self, Nx, Ny, Nz, T, rho):
		# Coordinates:
		x = np.linspace(0, self.boxsize - self.boxsize/Nx, Nx) + self.boxsize/Nx/2
		xx, yy, zz = np.meshgrid(x, x, x)
		Coord = np.zeros((Nx * Ny * Nz, 3))
		Coord[:,0] = xx.flatten()
		Coord[:,1] = yy.flatten()
		Coord[:,2] = zz.flatten()

		# Velocities:
		Vel = np.zeros((Nx * Ny * Nz, 3))

		# Mass:	
		#Mass = np.zeros(Nx * Ny * Nz) + (rho * self.boxsize**3 / (Nx * Ny * Nz))
		Mass = np.full(Nx*Ny*Nz, rho)

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
		Coord += self.boxsize/2	
	
		# Mass:
		#Mass = np.zeros(N) + (rho * (4/3) * np.pi * R**3 / N)
		Mass = np.full(N, rho)

		# Velocities:
		if J.all() != None:
			# Project positions onto the rotation axis (plane):
			a = Coord - self.boxsize/2
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

	def write_ic_file(self, fileName, partTypes, massTable=None):
		# Number of particle types:
		nPartTypes = 6
	
		with h5py.File(f'{self.path + fileName}.hdf5', 'w') as f:
			# Write particle types and their respective fields to hdf5-file:
			for ptype in partTypes.keys():
				g = f.create_group(ptype)
				for field in partTypes[ptype]:
					g[field] = partTypes[ptype][field]
		
			# Set particle counts:
			nPart = np.zeros(nPartTypes, dtype='int64')
			for ptype in partTypes.keys():
				pnum = int(ptype[-1])
				nPart[pnum] = len(partTypes[ptype]['ParticleIDs'])
		
			# Header:
			h = f.create_group('Header')
			h.attrs['BoxSize'] = self.boxsize
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

		return 0

	def check_ic(self, file):
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
	
			print('\nPartType1 keys:')
			for i in f['PartType1']:
				print(i)
	
			Coord_gas  = f['PartType0']['Coordinates'][:]
			Vel_gas    = f['PartType0']['Velocities'][:]
			Mass_gas   = f['PartType0']['Masses'][:]
			Coord_sink = f['PartType1']['Coordinates'][:]
			IntErg_gas = f['PartType0']['InternalEnergy'][:]

		# Plot:
		print('\nSimple plot of ICs:')
		fig, ax = plt.subplots(figsize=(6, 6))
		ax.scatter(Coord_gas[:,0], Coord_gas[:,1], marker='.', c='b')
		ax.scatter(Coord_sink[:,0], Coord_sink[:,1], marker='o', c='r')
		ax.set_aspect('equal')
		imgcat(fig)

		return 0

	def collapsing_sphere(self, Nx, Ny, Nz, N_sphere, T_grid, T_sphere, rho_grid, rho_sphere, R_sphere, J_sphere, Coord_sink, Vel_sink, M_sink, check=False):
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

		# Write IC file:
		self.write_ic_file(fileName='collapsing_sphere', partTypes={'PartType0' : ic_gas, 'PartType1' : ic_sink})
		print('Initial conditions generated!\n')

		# Check collapsing sphere properties:
		if check == True:
			print('Doing some quick esimates...')
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

			print('\nChecking generated HDF5-file...\n')
			self.check_ic(f'{self.path}/collapsing_sphere.hdf5')

		return 0


# -------------- Write IC-file
icetap = IceTap(path='/home/astro/jpeterss/ICGA/', boxsize=40, ulength=3.08567758e18, umass=1.9891e33, uvel=1e5)
icetap.collapsing_sphere(Nx=20, Ny=20, Nz=20, N_sphere=10000, T_grid=1.0e4, T_sphere=1.0e2, rho_grid=1.477e-2, rho_sphere=7.385, R_sphere=10, J_sphere=np.array([0, 0, 3.241e2]), 
					   Coord_sink=[20, 20, 20], Vel_sink=[0.0, 0.0, 0.0], M_sink=1e2, check=True)



'''
1 g/cm3 = 1.47705e22 Msol/pc3
1 cm2/s = 3.08568e23 pc km/s
Density: 1.0e-24 g/cm3 --> 1.477e-2
Density: 5.0e-22 g/cm3 --> 7.385
Specific angular momentum: 1e26 cm2/s --> 3.241e2
'''

# -------------- End of file


