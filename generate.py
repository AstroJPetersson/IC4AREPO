#!/home/astro/jpeterss/anaconda3/bin/python

# -------------- Packages
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import h5py
from imgcat import imgcat

from ic4arepo import write_ic_file, check_ic_file, SphereCollapse, GalaxyBarred


# -------------- Some good unit conversions
# 1 g/cm3 = 1.47705e22 Msol/pc3
# 1 cm2/s = 3.08568e23 pc km/s
# Density: 1.0e-24 g/cm3 --> 1.477e-2
# Density: 5.0e-22 g/cm3 --> 7.385
# Specific angular momentum: 1e26 cm2/s --> 3.241e2


# -------------- Generate IC-file
# Initialize argparse:
parser = argparse.ArgumentParser(description='ICE : Initial Conditions gEnerator for arepo', usage='icgenerate.py [options] model', 
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser._actions[0].help='Show this help message and exit'

# Positional arguments:
parser.add_argument('model', help='The model you would like to generate initial conditions for')

# Read arguments from the command line:
args = parser.parse_args()

# Generate IC-file of the given model:
if args.model == 'spherecollapse':
    C = SphereCollapse(boxSize=40, ulength=3.08567758e18, umass=1.9891e33, uvel=1e5)
    C.icgenerate(N_sphere=int(1e4), T_sphere=1e2, rho_sphere=7.385, R_sphere=10, J_sphere=np.array([0, 0, 3.241e2]), 
                 Coord_sink=[20, 20, 20], Vel_sink=[0, 0, 0], M_sink=1e2, 
                 N_grid=int(1e3), T_grid=1e4, rho_grid=1.477e-2,
                 filename='spherecollapse', savepath='/home/astro/jpeterss/IC4AREPO/ICs/', properties=True, check=True)

if args.model == 'galaxybarred':
    G = GalaxyBarred(boxSize=500, ulength=3.08567758e20, umass=1.9891e33, uvel=1e5)
    G.icgenerate(N=int(2e6), runs=2, wait='auto', jobDir='/hpcstorage/jpeterss/galaxybarred/meshrelax/', saveDir='/home/astro/jpeterss/IC4AREPO/ICs/')


# -------------- End of file


