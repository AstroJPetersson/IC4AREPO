#!/home/astro/jpeterss/anaconda3/bin/python

# -------------- Packages
import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from imgcat import imgcat

# ------------- IC modify class
class IceTapModify:
	'''
	ICETAPMOD : Initial Conditions gEneraTor for ArePo - MODifier
	'''
	def __init__(self, icfile):
		# Which file to modify:
		self.icfile = icfile

		# Constants:
		self.kb = 1.3807e-16     # [erg K^-1]
		self.G  = 6.6726e-8      # [dyne cm^2 g^-2]
		self.mp = 1.6726e-24     # [g]
		self.pc = 3.08567758e18  # [cm]

	def header(self, Time):
		with h5py.File(self.icfile, 'r+') as f:
			print(f['Header'].attrs['Time'])
			f['Header'].attrs['Time'] = Time
			

		return 0

	def check_ic(self):
		# Check generated IC-file:
		with h5py.File(self.icfile, 'r') as f:
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

# -------------- Modify IC-file
icma = ICMA(icfile='snap_002.hdf5')
icma.header(Time=0) 
icma.check_ic()

# -------------- End of file


