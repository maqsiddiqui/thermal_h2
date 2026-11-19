#!/bin/bash
#SBATCH --job-name=croft_h2_fv
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=100GB
#SBATCH --time=72:00:00
# Croft H2 storage, full grid, finite volume. Edit PFDIR and the input name; submit with: sbatch run_h2_hpc.sh
# To continue after the time limit: add --recover to the mpirun line and submit again.
source ~/set_env.sh
export PFDIR=$SCRATCH3DIR/moose/modules/porous_flow      # the porous_flow-opt built with the add-ons
cd "$SLURM_SUBMIT_DIR"
mpirun ${PFDIR}/porous_flow-opt -i croft_h2_fv.i > croft_h2_fv.out 2>&1
