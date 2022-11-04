import numpy as np
import h5py

class WriteHdf5:
	'''
	DOCSTRING
	'''
	def __init__(self, fileName, savePath, boxSize, partTypes, massTable=None):
		self.fileName  = fileName
		self.savePath  = savePath
		self.boxSize   = boxSize
		self.partTypes = partTypes
		self.massTable = massTable

	def write_ic_file(self):
		# Number of particle types:
		nPartTypes = 6
	
		with h5py.File(f'{self.savePath + self.fileName}.hdf5', 'w') as f:
			# Write particle types and their respective fields to hdf5-file:
			for ptype in self.partTypes.keys():
				group = f.create_group(ptype)
				for field in self.partTypes[ptype]:
					group[field] = self.partTypes[ptype][field]
		
			# Set particle counts:
			nPart = np.zeros(nPartTypes, dtype='int64')
			for ptype in self.partTypes.keys():
				pNum = int(ptype[-1])
				nPart[pNum] = len(self.partTypes[ptype]['ParticleIDs'])
		
			# Header:
			h = f.create_group('Header')
			h.attrs['BoxSize'] = self.boxSize
			h.attrs['NumFilesPerSnapshot'] = 1
			h.attrs['NumPart_ThisFile'] = nPart
			h.attrs['NumPart_Total'] = nPart
			h.attrs['NumPart_Total_HighWord'] = np.zeros(nPartTypes)
			
			for k in ['Time', 'Redshift', 'Omega0', 'OmegaLambda', 'HubbleParam']:
				h.attrs[k] = 0.0
			for k in ['Sfr', 'Cooling', 'StellarAge', 'Metals', 'Feedback']:
				h.attrs[f'Flag_{k}'] = 0
			h.attrs['Flag_DoublePrecision'] = 1

			if self.massTable is not None:
				h.attrs['MassTable'] = self.massTable
			else:
				h.attrs['MassTable'] = np.zeros(nPartTypes)

		print(f'Initial conditions file saved at: {self.savePath} as {self.fileName}.hdf5')

		return 0


