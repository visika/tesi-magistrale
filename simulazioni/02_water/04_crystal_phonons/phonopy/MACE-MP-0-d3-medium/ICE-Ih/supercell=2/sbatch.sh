#!/usr/bin/env sh
#SBATCH --job-name=phonopy
#SBATCH --partition=gpus
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --gpus=1
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

python3 ./phonons.py --model medium --d3 --system structs/ICE13/Ih/POSCAR --supercell 2 --device cuda --y_min -0.5 --grid 5 5 5

date
echo "FINISHED!!!"
rm $ALERT
