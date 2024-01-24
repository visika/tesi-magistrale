#!/usr/bin/bash
#SBATCH --job-name=Mollo
#SBATCH --nodes=1               # Number of nodes
#SBATCH --tasks-per-node=1     # Number of MPI ranks per node
##SBATCH --exclusive
##SBATCH --exclude=ibiscohpc-wn04,ibiscohpc-wn06,ibiscohpc-wn07
#SBATCH --cpus-per-task=1
##SBATCH --mem=1400000
#SBATCH --partition=sequential

# Requires: INCAR, INPHON, KPOINTS, POSCAR, POTCAR, SPOSCAR

### Parte iniziale dello script, con configurazioni iniziali. ###

ulimit -s unlimited

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

source /nfsexports/SOFTWARE/intel/oneapi/setvars.sh
export MPI="mpirun"
echo "     using $(which $MPI)"

# Prendi il percorso di VASP
export VASP="/lustre/home/tccourse/vasp46-da/vasp"

# Prendi il percorso di PHON
export PHON="/lustre/home/tccourse/Phon/src/phon"

echo "Il mio programma sta iniziando..."

# Passa alla cartella di esecuzione dello script
# Molto probabilmente è quella attuale
cd "$SLURM_SUBMIT_DIR" || exit

### Da qui inizia lo script con i comandi specifici. ###

$PHON >phon.out
