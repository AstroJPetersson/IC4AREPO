import numpy as np
import matplotlib.pyplot as plt
import h5py
from imgcat import imgcat

def check_ic_file(file):
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

		print(f'\nMass of PartType0 particles [msol]: {np.unique(Mass_gas)}')

	# Plot:
	print('\nSimple plot of ICs:')
	fig, ax = plt.subplots(figsize=(6, 6))
	ax.scatter(Coord_gas[:,0], Coord_gas[:,1], marker='.', c='b')
	ax.scatter(Coord_sink[:,0], Coord_sink[:,1], marker='o', c='r')
	ax.set_aspect('equal')
	imgcat(fig)

	return 0


