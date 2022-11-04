import numpy as np
import matplotlib.pyplot as plt
import h5py
from imgcat import imgcat

class Galaxy:
	'''
	DOCSTRING
	'''
	def __init__(self, boxSize, ulength, umass, uvel):
		# Constants:
		self.kb   = 1.3807e-16     # [erg K^-1]
		self.G    = 6.6726e-8      # [dyne cm^2 g^-2]
		self.mp   = 1.6726e-24     # [g]
		self.pc   = 3.08567758e18  # [cm]
		self.Msol = 1.9891e33      # [g] 
		self.kpc  = 3.08567758e21  # [cm]


	def galaxy(self):
		'''
		Galaxy model.
		'''

		return 0


	def mesh_relax(self):
		'''
		AREPO mesh relax. 
		'''

		return 0

	def density_profile(self, pos):
		'''
		Density profile for the galaxy model (source: Robin Tress). 
		'''
		Sigma0 = 50  * self.Msol/self.pc/self.pc / self.ucoldens
		zd     = 85  * self.pc                   / self.ulength
		Rm     = 1.5 * self.kpc                  / self.ulength
		Rd     = 7   * self.kpc                  / self.ulength	
		R_cut  = 10  * self.kpc                  / self.ulength

		x, y, z = pos.T - self.boxSize/2
		R = np.sqrt(x**2 + y**2)
		density = Sigma0/(4*zd) * np.exp(-Rm/R - R/Rd) / np.cosh(z/(2*zd))**2
		density[R>R_cut] = 1e-28 / self.udens

		return density


