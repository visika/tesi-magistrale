#!/usr/bin/env bash
#SBATCH --job-name=vib
#SBATCH --nodes=1             # Number of nodes
#SBATCH --tasks-per-node=4    # Number of MPI ranks per node
#SBATCH --cpus-per-task=1
#SBATCH --partition=parallel
#SBATCH --propagate=NONE

ulimit -s unlimited

# export OMP_NUM_THREADS=4
# export MKL_NUM_THREADS=1

# source /nfsexports/SOFTWARE/intel/oneapi/setvars.sh
# export MPI="mpirun"
# echo "     using $(which $MPI)"

# Prendi il percorso di VASP
# export passa l'assegnazione anche ai processi figli
# export VASP="/lustre/home/tccourse/vasp46-da/vasp"
# export VASPGAMMA="/lustre/home/tccourse/vasp46-da.gamma/vasp"
# Prendi il percorso di PHON
# export PHON="/lustre/home/tccourse/Phon/src/phon"
# export RUNPHON="$SLURM_SUBMIT_DIR/runphon"

# Passa alla cartella di esecuzione dello script
cd "$SLURM_SUBMIT_DIR" || exit

### Da qui inizia lo script con i comandi specifici. ###
ALERT=000-THIS_IS_RUNNING
touch $ALERT

source /lustre/home/mmollo/setupconda.sh
conda activate mmollo-lammps-env

export VASP_PP_PATH=/lustre/home/mmollo/vasp_pp

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/nvidia/hpc_sdk/Linux_x86_64/22.3/comm_libs/openmpi4/openmpi-4.0.5/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/extras/qd/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin

echo Using "$(python -V)"
python test_al.py

echo "FINISHED!!!"
rm $ALERT
