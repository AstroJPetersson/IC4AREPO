#!/home/astro/jpeterss/anaconda3/bin/python

# -------------- Packages
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py
from imgcat import imgcat

from src import CollapsingSphere, WriteHdf5, check_icfile

# -------------- Generate IC-file
c = CollapsingSphere(boxSize=40, ulength=3.08567758e18, umass=1.9891e33, uvel=1e5)
ic_gas, ic_sink = c.icgenerate(Nx=20, Ny=20, Nz=20, N_sphere=10000, T_grid=1.0e4, T_sphere=1.0e2, rho_grid=1.477e-2, 
							   rho_sphere=7.385, R_sphere=10, J_sphere=np.array([0, 0, 3.241e2]),  
							   Coord_sink=[20, 20, 20], Vel_sink=[0.0, 0.0, 0.0], M_sink=1e2)

w = WriteHdf5(boxSize=40, fileName='collapsing_sphere', savePath='/home/astro/jpeterss/ICETAP/ICs/', 
			  partTypes={'PartType0' : ic_gas, 'PartType1' : ic_sink})
w.write_ic_file()

check_icfile(icfile='ICs/collapsing_sphere.hdf5')

# 1 g/cm3 = 1.47705e22 Msol/pc3
# 1 cm2/s = 3.08568e23 pc km/s
# Density: 1.0e-24 g/cm3 --> 1.477e-2
# Density: 5.0e-22 g/cm3 --> 7.385
# Specific angular momentum: 1e26 cm2/s --> 3.241e2

# -------------- End of file


