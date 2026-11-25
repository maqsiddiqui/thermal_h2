#!/bin/bash
#SBATCH --job-name=croft_hist_fv
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=100GB
#SBATCH --time=72:00:00
# Croft history match, full grid, finite volume (hysteresis + fault connections).
# Edit EXE and INPUT; submit with: sbatch run_hist_hpc.sh
# To continue after the time limit: set RECOVER=--recover below and submit again.
source ~/set_env.sh
EXE=/scratch3/sid036/thermal_h2/thermal_h2-opt     # your MOOSE app, built with the add-ons
INPUT=croft_full_fv_nnc
RECOVER=
cd "$SLURM_SUBMIT_DIR"
mpirun ${EXE} -i ${INPUT}.i ${RECOVER} >> ${INPUT}.out 2>&1
