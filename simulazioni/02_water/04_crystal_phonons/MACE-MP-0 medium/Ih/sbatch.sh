#!/usr/bin/env sh
#SBATCH --job-name=phonIh
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=8
#SBATCH --propagate=NONE

echo Starging at $(date)

ulimit -s unlimited

source /nfsexports/SOFTWARE/intel/oneapi/setvars.sh
export MPI="mpirun"
echo "     using $(which $MPI)"

cd "$SLURM_SUBMIT_DIR" || exit

ALERT=000-THIS_IS_RUNNING
touch $ALERT

source $(poetry env info --path)/bin/activate

echo Using "$(python -V)"
echo Python path: "$(which python)"

python script.py

date
echo "FINISHED!!!"
rm $ALERT
