import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import time
import subprocess
from imgcat import imgcat

class Galaxy:
	'''
	DOCSTRING
	'''
	def __init__(self, boxSize, ulength, umass, uvel):
		# Simulation set-up:
		self.boxSize  = boxSize
		self.ulength  = ulength
		self.uvel     = uvel
		self.utime    = ulength/uvel
		self.udens    = umass/(ulength**3)
		self.ucoldens = umass/(ulength**2)
		self.uinterg  = uvel**2
		
		# Constants:
		self.kb   = 1.3807e-16     # [erg K^-1]
		self.G    = 6.6726e-8      # [dyne cm^2 g^-2]
		self.mp   = 1.6726e-24     # [g]
		self.pc   = 3.08567758e18  # [cm]
		self.Msol = 1.9891e33      # [g] 
		self.kpc  = 3.08567758e21  # [cm]

	def galaxymodel(self, N):
		'''
		Artificial ICs for a Milky Way-like galaxy. 
		'''
		# Positions:
		pos = self.boxSize * np.random.rand(N, 3)

		# Velocities: 
		# (UNDER CONSTRUCTION)
		vel = np.zeros((N, 3))

		# Densities:
		dens = self.density_profile(pos)
		
		# Internal energy:
		interg = np.zeros(N)

		# Particle IDs:
		iord = np.arange(1, N+1)

		# Initial conditions:
		ic_galaxy = {
			'Coordinates'    : pos,
			'Velocities'     : vel,
			'Masses'         : dens,
			'InternalEnergy' : interg,
			'ParticleIDs'    : iord
		}

		return ic_galaxy

	def mesh_relax(self, jobDir, wait='manuel', sleep=120):
		initDir = os.getcwd()
		os.chdir(jobDir)

		print('Submitting and waiting for job to be completed...')
		os.system('sbatch job.sh')
		
		# 'sleep' wait method:
		if wait == 'sleep':
			if os.path.isfile('Arepo.out'):
				with open('Arepo.out') as A:
					a = A.readlines()
					line_cut = a[-4][13:-10]
					last_run = float(line_cut)
					sleep = 1.3*last_run
					print(f'Last code run: {last_run} s | Estimated code run: {sleep} s')
			time.sleep(sleep)
		# 'auto' wait method:
		elif wait == 'auto':
			print('(you are running with wait=\'auto\': the code will now automatically check whether the final snapshot has been produced, and not proceed until it has)')
			while not os.path.exists(f'{jobDir}/output/snap_001.hdf5'):
				time.sleep(1)
		# 'manual' wait method:
		else:
			print('(please keep track of the submitted job manually)')
			proc = 'n'
			while proc == 'n':
				proc = input('Is the job completed? (y/n): ')
		
		os.chdir(initDir)

		return 0

	def icmodify(self, jobDir):
		# Move final snapshot from previous mesh-relax and rename it to galaxy.hdf5
		os.system(f'mv {jobDir}/output/snap_001.hdf5 {jobDir}/galaxy.hdf5')

		# Modify the new IC-file with the density profile:
		with h5py.File(f'{jobDir}/galaxy.hdf5', 'r+') as f:
			pos = f['PartType0']['Coordinates'][:]
			dens_upd = self.density_profile(pos)
			f['PartType0']['Masses'][:] = dens_upd

		return 0

	def icgenerate(self, N, runs, jobDir, saveDir, wait='manual', sleep=120):
		'''
		Generate IC-file using AREPO mesh-relax. 
		''' 	
		from .write import write_ic_file

		# Generate galaxy model and wrtie ic-file:
		print('Generating artifical ICs for a Milky Way-like galaxy...')
		ic_galaxy = self.galaxymodel(N)
		write_ic_file(name='galaxy', path=jobDir, boxSize=self.boxSize, partTypes={'PartType0' : ic_galaxy})
		print('DONE!\n')

		# Run and repeat mesh relaxation on ic-file:
		print(f'Initialising mesh-relaxtion on the generated IC-file:')
		for i in range(1, runs+1):
			print(f'Mesh-relaxation {i}/{runs}: Running...')
			self.mesh_relax(jobDir=jobDir, wait=wait, sleep=sleep)
			print(f'Mesh-relaxation {i}/{runs}: DONE!\n')

			if i < runs:
				print('Preparing IC-file for next mesh-relax...')
				self.icmodify(jobDir=jobDir)
				print('Preparation done!\n')
			elif i == runs:
				print(f'All interations of mesh-relaxation done! New IC-file can be found at: \'{saveDir}\' as \'galaxy.hdf5\'')
				os.system(f'mv {jobDir}/output/snap_001.hdf5 {saveDir}/galaxy.hdf5')

		return 0

	def density_profile(self, pos):
		'''
		Density profile for the galaxy model (ref: Robin Tress). 
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


