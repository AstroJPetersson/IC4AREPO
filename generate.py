#!/home/astro/jpeterss/anaconda3/bin/python

# -------------- Packages:
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import h5py
from imgcat import imgcat

from src import CloudCollapse, GalaxyBarred, Streams


# -------------- Units:
pc   = 3.08567758e18  #[cm]
kpc  = 3.08567758e21  #[cm]
Msol = 1.9891e33      #[g]
kmps = 1e5            #[cm/s]
Myr  = 1e6 * 365.25 * 24 * 60 * 60  #[s]

# -------------- Arguments:
# Initialize argparse:
parser = argparse.ArgumentParser(description='IC4AREPO : Initial Conditions 4 AREPO', 
                                 usage='icgenerate.py [options] model', 
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser._actions[0].help='Show this help message and exit'

# Positional arguments:
parser.add_argument('model', help='Model selection')

# Read arguments from the command line:
args = parser.parse_args()


# -------------- Generate IC-file
print('\nWelcome to IC4AREPO:')
print(f'  * Model: \'{args.model}\'')
print('  * Starting to generate initial conditions... please wait...')

# -------------- Uniform Cloud
if args.model == 'uniformcloud':
    # Parameters:
    N         = 10000000
    T         = 1e2
    rho       = 5e-22
    R         = 10 * pc 
    P         = 0.0001 * Myr
    rot       = np.array([0, 0, 1])
    pos_cloud = np.array([50, 50, 50]) * pc
    vel_cloud = np.array([0, 0, 0]) * kmps

    pos_sink  = np.array([50, 50, 50]) * pc
    vel_sink  = np.array([0, 0, 0])
    mass_sink = 1e7 * Msol
    
    # Generate ICs:
    C = CloudCollapse(boxSize=100 * pc, ulength=3.08567758e18, umass=1.9891e33, uvel=1e5)
    C.icgenerate_uniform_cloud(N=N, pos_cloud=pos_cloud, vel_cloud=vel_cloud, R=R, rho=rho, T=T, P=P, 
                               rot=rot, pos_sink=pos_sink, vel_sink=vel_sink, mass_sink=mass_sink, 
                               filename='cloudcollapse', savepath='/home/astro/jpeterss/IC4AREPO/ICs/', 
                               check=True, relax=True, N_relax=5, wait='auto')

# -------------- R2 Cloud
if args.model == 'r2cloud':
    # Parameters:
    N         = 10000000
    T         = 1e2
    rho       = 5e-22
    R         = 10 * pc 
    pos_cloud = np.array([50, 50, 50]) * pc
    vel_cloud = np.array([0, 0, 0]) * kmps

    pos_sink  = np.array([50, 50, 50]) * pc
    vel_sink  = np.array([0, 0, 0])
    mass_sink = 1e7 * Msol
    
    # Generate ICs:
    C = CloudCollapse(boxSize=100 * pc, ulength=3.08567758e18, umass=1.9891e33, uvel=1e5)
    C.icgenerate_r2_cloud(N=N, pos_cloud=pos_cloud, vel_cloud=vel_cloud, R=R, rho=rho, T=T, 
                          pos_sink=pos_sink, vel_sink=vel_sink, mass_sink=mass_sink, 
                          filename='cloudcollapse', savepath='/home/astro/jpeterss/IC4AREPO/ICs/', 
                          check=True, relax=True, N_relax=5, wait='auto')

# -------------- Two Uniform Clouds
if args.model == 'twouniformclouds':
    # Parameters:
    N          = 200000
    T1, T2     = 1e2, 1e2
    rho1, rho2 = 1e-21, 1e-21
    R1, R2     = 10 * pc, 10 * pc 
    pos_cloud1, pos_cloud2 = np.array([20, 20, 50]) * pc, np.array([80, 80, 50]) * pc
    vel_cloud1, vel_cloud2 = np.array([0, 0, 0]), np.array([0, 0, 0]) 

    pos_sink  = np.array([50, 50, 50]) * pc
    vel_sink  = np.array([0, 0, 0])
    mass_sink = 1e7 * Msol
    
    # Generate ICs:
    C = CloudCollapse(boxSize=100 * pc, ulength=3.08567758e18, umass=1.9891e33, uvel=1e5)
    C.icgenerate_twoclouds(N=N, pos_cloud1=pos_cloud1, pos_cloud2=pos_cloud2, vel_cloud1=vel_cloud1, 
                           vel_cloud2=vel_cloud2, R1=R1, R2=R2, rho1=rho1, rho2=rho2, T1=T1, T2=T2, 
                           pos_sink=pos_sink, vel_sink=vel_sink, mass_sink=mass_sink, 
                           filename='cloudcollapse', savepath='/home/astro/jpeterss/IC4AREPO/ICs/', 
                           check=True, relax=True, N_relax=1, wait='manual')

# -------------- Two Opposite Streams
if args.model == 'streams':
    # Parameters:
    
    S = Streams

# -------------- Barred Galaxy
if args.model == 'galaxybarred':
    G = GalaxyBarred(boxSize=500, ulength=3.08567758e20, umass=1.9891e33, uvel=1e5)

# -------------- End of file
print('  * The End, happy simulation :)\n')


