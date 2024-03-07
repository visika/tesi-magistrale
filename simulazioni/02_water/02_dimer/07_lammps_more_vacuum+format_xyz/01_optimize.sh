#!/usr/bin/env bash
#SBATCH --job-name=dimer
#SBATCH --nodes=1            # Number of nodes
#SBATCH --ntasks-per-node=1  # Number of MPI ranks per node
#SBATCH --cpus-per-task=4    # = 1 in sequential partition
#SBATCH --partition=parallel

ulimit -s unlimited

# export OMP_NUM_THREADS=1 # Corrisponde a --cpus-per-task sopra
# export MKL_NUM_THREADS=1

# source /nfsexports/SOFTWARE/intel/oneapi/setvars.sh
# export MPI="mpirun"
# echo "     using $(which $MPI)"

# Prendi il percorso di VASP
# export passa l'assegnazione anche ai processi figli
export VASP="/lustre/home/tccourse/vasp46-da/vasp"
export VASPGAMMA="/lustre/home/tccourse/vasp46-da.gamma/vasp"
# Prendi il percorso di PHON
export PHON="/lustre/home/tccourse/Phon/src/phon"
export RUNPHON="$SLURM_SUBMIT_DIR/runphon"

# Passa alla cartella di esecuzione dello script
cd "$SLURM_SUBMIT_DIR" || exit

### Da qui inizia lo script con i comandi specifici. ###
ALERT=000-THIS_IS_RUNNING
touch $ALERT

source /lustre/home/mmollo/setupconda.sh
conda activate mmollo-lammps-env

echo Using "$(python -V)"
python 01_optimize.py

echo "FINISHED!!!"
rm $ALERT
