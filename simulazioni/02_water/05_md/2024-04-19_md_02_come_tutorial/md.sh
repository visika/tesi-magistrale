#!/usr/bin/env bash
#SBATCH --job-name=MD
#SBATCH --partition=gpus
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --gpus=1
#SBATCH --propagate=NONE

ulimit -s unlimited

source /nfsexports/SOFTWARE/intel/oneapi/setvars.sh
export MPI="mpirun"
echo "     using $(which $MPI)"

# Questo serve per poetry
export PATH=$PATH:/lustre/home/mmollo/.local/bin

cd "$SLURM_SUBMIT_DIR" || exit

ALERT=000-THIS_IS_RUNNING
touch $ALERT

source "$(poetry env info --path)/bin/activate"

echo Using "$(python -V)"
echo Python path: "$(which python)"

python md.py

echo "FINISHED!!!"
rm $ALERT
