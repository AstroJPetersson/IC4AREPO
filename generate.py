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
# 1 cm2/s = 3.24078e-24 pc km/s

# Density: 1.0e-24 g/cm3 --> 1.477e-2
# Density: 5.0e-22 g/cm3 --> 7.385
# Specific angular momentum: 1e26 cm2/s --> 3.24e2

pc   = 3.08567758e18  #[cm]
kpc  = 3.08567758e21  #[cm]
Msol = 1.9891e33      #[g]
kmps = 1e5            #[cm/s]

# -------------- Generate IC-file
# Initialize argparse:
parser = argparse.ArgumentParser(description='IC4AREPO : Initial Conditions 4 AREPO', usage='icgenerate.py [options] model', 
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser._actions[0].help='Show this help message and exit'

# Positional arguments:
parser.add_argument('model', help='Model for which you will generate initial conditions of')

# Read arguments from the command line:
args = parser.parse_args()

# Generate IC-file:
print('\nWelcome to IC4AREPO:')
print(f'  * Model: {args.model}')
print('  * Starting to generate initial conditions...')

if args.model == 'spherecollapse':
    C = SphereCollapse(boxSize=40 * pc, ulength=3.08567758e18, umass=1.9891e33, uvel=1e5)
    
    N_sphere, N_grid     = 100000, 10000
    T_sphere, T_grid     = 1e2, 1e4
    rho_sphere, rho_grid = 1e-21, 1e-24
    R_sphere             = 10 * pc 
    J_sphere             = np.array([0, 0, 1e27])
    
    pos_sink  = np.array([20, 20, 20]) * pc
    vel_sink  = np.array([0, 0, 0])
    mass_sink = 1e7 * Msol
    
    C.icgenerate(N_sphere=N_sphere, T_sphere=T_sphere, rho_sphere=rho_sphere, R_sphere=R_sphere, J_sphere=J_sphere, 
                 pos_sink=pos_sink, vel_sink=vel_sink, mass_sink=mass_sink, N_grid=N_grid, T_grid=T_grid, rho_grid=rho_grid,
                 filename='spherecollapse', savepath='/home/astro/jpeterss/IC4AREPO/ICs/', jeans_check=True, check=True)

if args.model == 'galaxybarred':
    G = GalaxyBarred(boxSize=500, ulength=3.08567758e20, umass=1.9891e33, uvel=1e5)
    G.debug()
    #G.icgenerate(N_gal=int(1e6), N_grid=int(1e6), runs=5, wait='auto', jobDir='/hpcstorage/jpeterss/galaxybarred/meshrelax/', 
    #             saveDir='/home/astro/jpeterss/IC4AREPO/ICs/', check=True)


print('  * The End, happy simulation :)\n')

# -------------- End of file


