#!/home/astro/jpeterss/anaconda3/bin/python

# -------------- packages
import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from imgcat import imgcat

# -------------- units
pc      = 3.08567758e18
msun    = 1.9891e33

# -------------- internal units
ulength  = pc
umass    = msun
uvel     = 1e5
udensity = umass / ulength**3
uenergy  = uvel**2
utime    = ulength / uvel
uangmom  = ulength * uvel

# -------------- constants
kb = 1.3807e-16               # [erg K^-1]
G  = 6.6726e-8                # [dyne cm^2 g^-2]
mH = 1.00794 * 1.6605402e-24  # [g]
mp = 1.6726e-24               # [g]

# -------------- simulation set-up
simdir = '/home/astro/jpeterss/ICs'
icfile = simdir + '/icsphere.hdf5'

# Grid:
boxSize    = 40 * pc / ulength
Nx, Ny, Nz = 20, 20, 20
T0Grid     = 1.0e4
rho0Grid   = 1.0e-24 / udensity

# Sphere:
rSphere    = 10 * pc / ulength
NSphere    = 10000
T0Sphere   = 1.0e2
rho0Sphere = 5.0e-22 / udensity
JSphere    = np.array([0, 0, 1e26]) / uangmom

# Sink:
posSink = np.array([20.0, 20.0, 20.0])
velSink = np.array([0.0, 0.0, 0.0])
mSink   = 1e2

# -------------- some estimates
# Sound speed:
gamma = 5/3
cs    = ((gamma * kb * T0Sphere) / (mH))**(1/2)

# Jeans length:
jl = cs / (G * rho0Sphere * udensity)**(1/2)

# Free-fall time:
tff = ((3*np.pi)/(32*G*rho0Sphere*udensity))**(1/2)

print(f'Jeans length [pc] : {jl/pc}')
print(f'Sound speed [km s^-1] : {cs/(1e5)}')
print(f'Free-fall time [Myr]: {tff/(1e6 * 365.25 * 24 * 60 * 60)}')

# -------------- functions
def generate_grid(Nx, Ny, Nz, boxSize, rho0, T0):
	pos = np.zeros((Nx * Ny * Nz, 3))
	vel = np.zeros((Nx * Ny * Nz, 3))
	
	# Coordinates:
	x = np.linspace(0, boxSize - boxSize/Nx, Nx) + boxSize/Nx/2
	y = np.linspace(0, boxSize - boxSize/Ny, Ny) + boxSize/Ny/2
	z = np.linspace(0, boxSize - boxSize/Nz, Nz) + boxSize/Nz/2
	
	xx, yy, zz = np.meshgrid(x, y, z)
	
	pos[:,0] = xx.flatten()
	pos[:,1] = yy.flatten()
	pos[:,2] = zz.flatten()

	# Mass:	
	mass = np.zeros(Nx * Ny * Nz) + (rho0 * boxSize**3 / (Nx * Ny * Nz))

	# Internal Energy:
	mu = 1
	u = np.full(Nx * Ny * Nz, (3 * kb * T0) / (2 * mu * mH) / uenergy)
	
	# ICs:
	ic_grid = {
		'Coordinates': pos,
		'Velocities': vel,
		'Masses': mass,
		'InternalEnergy': u,
	}
	
	return ic_grid

def generate_sphere(N, boxSize, R, rho0, T0, J=None):
	# Coordinates:
	pos = np.zeros(shape=(N, 3))
	mask = np.full(N, False)
	while len(pos[mask]) < N:
		pos[:,0][~mask] = 2*R*np.random.random(N - len(pos[mask])) - R
		pos[:,1][~mask] = 2*R*np.random.random(N - len(pos[mask])) - R
		pos[:,2][~mask] = 2*R*np.random.random(N - len(pos[mask])) - R
		mask = (np.linalg.norm(pos, axis=1) >= 0.1) * (np.linalg.norm(pos, axis=1) <= R)
	pos += boxSize/2	
	
	# Mass:
	mass = np.zeros(N) + (rho0 * (4/3) * np.pi * R**3 / N)
	
	 # Velocities:
	if J.all() != None:
		# Project positions onto the rotation axis (plane):
		a = pos - boxSize/2
		a1 = (np.dot(a, J) / (np.linalg.norm(J)**2))[:,np.newaxis] * J
		a2 = a - a1
		# Calculate velocities given the total (specific) angular momentum:
		v = (np.linalg.norm(J)/N) * np.linalg.norm(a2, axis=1)
		perp_vec = np.cross(J, a2)
		vel = v[:,np.newaxis] * perp_vec / (np.linalg.norm(perp_vec, axis=1)[:,np.newaxis])
	else:
		vel = np.zeros((N, 3))
	
	# Internal Energy:
	mu = 1
	u = np.full(N, (3 * kb * T0) / (2 * mu * mH) / uenergy)
	
	# ICs:
	ic_sphere = {
		'Coordinates': pos,
		'Velocities': vel,
		'Masses': mass,
		'InternalEnergy': u,
	}
	    
	return ic_sphere

def generate_sink(x, y, z, vx, vy, vz, m):
	ic_sink = {
		'Coordinates' : np.array([[x, y, z]]),
		'Velocities' : np.array([[vx, vy, vz]]),
		'Masses' : np.array([m]),
	}
	
	return ic_sink

def write_ic_file(fileName, partTypes, boxSize, massTable=None):
	# Number of particle types:
	nPartTypes = 6
	
	with h5py.File(fileName, 'w') as f:
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
		
		# Create header:
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
	
	return 0

# -------------- generate IC file
# Generate and combine initial conditions for PartType0:
ic_grid  = generate_grid(Nx=Nx, Ny=Ny, Nz=Nz, boxSize=boxSize, rho0=rho0Grid, T0=T0Grid)
ic_sphere = generate_sphere(N=NSphere, boxSize=boxSize, R=rSphere, rho0=rho0Sphere, T0=T0Sphere, J=JSphere)

ic_gas = {}

for k in ic_grid.keys():
	key_grid = ic_grid[k]
	key_sphere = ic_sphere[k]
	key = np.append(key_grid, key_sphere, axis=0)
	ic_gas[k] = key

# Generate initial conditions for PartType1:
ic_sink = generate_sink(x=posSink[0], y=posSink[1], z=posSink[2], vx=velSink[0], vy=velSink[1], vz=velSink[2], m=mSink)

# Generate ParticleIDs:
ic_gas['ParticleIDs']  = np.arange(0, NSphere + Nx * Ny * Nz) + 1
ic_sink['ParticleIDs'] = np.array([ic_gas['ParticleIDs'][-1] + 1])

# Write IC file:
write_ic_file(fileName=icfile, partTypes={'PartType0' : ic_gas, 'PartType1' : ic_sink}, boxSize=boxSize)
print('Initial conditions generated!')

# -------------- check IC file
with h5py.File(icfile, 'r') as f:
	print('HDF5 file keys:')
	for i in f.keys():
		print(i)
	
	print('\nHDF5 Header attributes:')
	for i in f['Header'].attrs:
		v = f['Header'].attrs[i]
		print(f'{i} : {v}')
	
	print('\nPartType0 keys:')
	for i in f['PartType0']:
		print(i)
	
	print('\nPartType1 keys:')
	for i in f['PartType1']:
		print(i)
	
	Coord_gas = f['PartType0']['Coordinates'][:]
	Vel_gas = f['PartType0']['Velocities'][:]
	Mass_gas = f['PartType0']['Masses'][:]
	Coord_sink = f['PartType1']['Coordinates'][:]
	IE_gas = f['PartType0']['InternalEnergy'][:]

# Plot ICs:	
fig, ax = plt.subplots(figsize=(6, 6))

ax.scatter(Coord_gas[:,0], Coord_gas[:,1], marker='.', c='b')
ax.scatter(Coord_sink[:,0], Coord_sink[:,1], marker='o', c='r')
ax.set_aspect('equal')
imgcat(fig)

# ------------- end of file


