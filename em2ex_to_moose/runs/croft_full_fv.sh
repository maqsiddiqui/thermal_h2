#!/bin/bash -l
#SBATCH --job-name=croft_full_fv
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=72:00:00
#SBATCH --mem=100GB
#SBATCH --mail-user=Mohammed.Siddiqui.csiro.au
#SBATCH --mail-type=ALL
#SBATCH --account=OD-246275

source ~/set_env.sh
cd /scratch3/sid036/problems/croft_fv
mpirun /scratch3/sid036/thermal_h2/thermal_h2-opt -i croft_full_fv.i > croft_full_fv.out
