#!/bin/bash
#SBATCH --job-name=croft_h2_fv
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=100GB
#SBATCH --time=72:00:00
# Croft H2 storage, full grid, finite volume. Edit EXE and the input name; submit with: sbatch run_h2_hpc.sh
# To continue after the time limit: add --recover to the mpirun line and submit again.
source ~/set_env.sh
EXE=/scratch3/sid036/thermal_h2/thermal_h2-opt     # your MOOSE app, built with the add-ons
cd "$SLURM_SUBMIT_DIR"
mpirun ${EXE} -i croft_h2_fv.i > croft_h2_fv.out 2>&1
