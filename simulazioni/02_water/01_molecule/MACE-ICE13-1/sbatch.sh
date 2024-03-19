#!/usr/bin/env sh
#SBATCH --job-name=H2O
#SBATCH --partition=gpus
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --gpus=1
#SBATCH --propagate=NONE

ulimit -s unlimited

source /nfsexports/SOFTWARE/intel/oneapi/setvars.sh
export MPI="mpirun"
echo "     using $(which $MPI)"

# Passa alla cartella di esecuzione dello script
cd "$SLURM_SUBMIT_DIR" || exit

### Da qui inizia lo script con i comandi specifici. ###
ALERT=000-THIS_IS_RUNNING
touch $ALERT

source $(poetry env info --path)/bin/activate

echo Using "$(python -V)"
echo Python path: "$(which python)"

python script.py

echo "FINISHED!!!"
rm $ALERT
