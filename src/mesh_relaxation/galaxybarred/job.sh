#!/bin/sh
#SBATCH --job-name relax
#SBATCH --error err.job%j
#SBATCH --output out.job%j
#SBATCH --nodes 1 
#SBATCH --exclusive
#SBATCH --partition p4
#SBATCH --time 0-24:00:00

module purge

. /usr/local/etc/setup-modules.sh

ml load foss/2021b
ml load GSL
ml load GMP
ml load Eigen
ml load HDF5

SIMDIR=$(pwd)

if [[ ! -d $SIMDIR/OUTPUT ]]; then
	mkdir OUTPUT
fi

mpiexec ./Arepo param.txt > Arepo.out


