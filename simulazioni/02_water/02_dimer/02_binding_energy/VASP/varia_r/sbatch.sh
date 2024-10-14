#!/usr/bin/env bash
#SBATCH --job-name=bind
#SBATCH --partition=gpus
#SBATCH --gpus=4
#SBATCH --ntasks-per-node=4
#SBATCH --propagate=NONE

ulimit -s unlimited

# Importa pacchetti per l'esecuzione di VASP
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/nvidia/hpc_sdk/Linux_x86_64/22.3/comm_libs/openmpi4/openmpi-4.0.5/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/extras/qd/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin
export vasp_std=/ibiscostorage/VirtualMatterLab/vasp_6.3.0_bin/vasp_std

# Imposta le variabili per l'uso di VASP con ASE
export ASE_VASP_COMMAND="mpirun -np 4 $vasp_std"
export VASP_PP_PATH=/lustre/home/mmollo/vasp_pp

cd "$SLURM_SUBMIT_DIR" || exit

ALERT=000-THIS_IS_RUNNING
touch $ALERT

# Questo serve per poetry
export PATH=$PATH:/lustre/home/mmollo/.local/bin
source "$(poetry env info --path)/bin/activate"

echo Using "$(python -V)"
echo Python path: "$(which python)"
python binding_energy.py

echo "FINISHED!!!"
rm $ALERT
