#!/usr/bin/env bash
#SBATCH --job-name=vib
#SBATCH --nodes=1             # Number of nodes
#SBATCH --tasks-per-node=1    # Number of MPI ranks per node
#SBATCH --cpus-per-task=1
#SBATCH --partition=sequential

ulimit -s unlimited

# Importa pacchetti per l'esecuzione di VASP
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/nvidia/hpc_sdk/Linux_x86_64/22.3/comm_libs/openmpi4/openmpi-4.0.5/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/extras/qd/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin
export vasp_std=/ibiscostorage/VirtualMatterLab/vasp_6.3.0_bin/vasp_std

# Imposta le variabili per l'uso di VASP con ASE
export ASE_VASP_COMMAND="mpirun $vasp_std"
export VASP_PP_PATH=/lustre/home/mmollo/vasp_pp

# Passa alla cartella di esecuzione dello script
cd "$SLURM_SUBMIT_DIR" || exit

### Da qui inizia lo script con i comandi specifici. ###
ALERT=000-THIS_IS_RUNNING
touch $ALERT

# Questo serve per poetry
export PATH=$PATH:/lustre/home/mmollo/.local/bin
source "$(poetry env info --path)/bin/activate"

echo Using "$(python -V)"
echo Python path: "$(which python)"

python 01_optimize.py

echo "FINISHED!!!"
rm $ALERT
