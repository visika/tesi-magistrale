#!/usr/bin/env bash
#SBATCH --job-name=Mollo
#SBATCH --nodes=1               # Number of nodes
#SBATCH --tasks-per-node=4     # Number of MPI ranks per node
##SBATCH --exclusive
##SBATCH --exclude=ibiscohpc-wn04,ibiscohpc-wn06,ibiscohpc-wn07
#SBATCH --cpus-per-task=1
##SBATCH --mem=1400000
#SBATCH --partition=parallel

### Parte iniziale dello script, con configurazioni iniziali. ###

ulimit -s unlimited

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

source /nfsexports/SOFTWARE/intel/oneapi/setvars.sh
export MPI="mpirun"
echo "     using $(which $MPI)"

# Prendi il percorso di VASP
# export passa l'assegnazione anche ai processi figli
export VASP="/lustre/home/tccourse/vasp46-da/vasp"
export VASPGAMMA="/lustre/home/tccourse/vasp46-da.gamma/vasp"

# Prendi il percorso di PHON
export PHON="/lustre/home/tccourse/Phon/src/phon"

export RUNPHON="$SLURM_SUBMIT_DIR/runphon"

echo "Il mio programma sta iniziando..."

# Passa alla cartella di esecuzione dello script
# Molto probabilmente è quella attuale
cd "$SLURM_SUBMIT_DIR" || exit

### Da qui inizia lo script con i comandi specifici. ###

ALERT=000-THIS_IS_RUNNING
touch $ALERT
./simultaneous.sh
echo "    FINISHED!!!"
rm $ALERT
