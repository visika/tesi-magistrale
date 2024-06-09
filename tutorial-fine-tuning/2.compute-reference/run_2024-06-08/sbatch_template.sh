#!/usr/bin/bash
#SBATCH --job-name=VASP
#SBATCH --partition=parallel
#SBATCH --ntasks-per-node=4 # Number of MPI ranks per node
#SBATCH --propagate=NONE

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/nvidia/hpc_sdk/Linux_x86_64/22.3/comm_libs/openmpi4/openmpi-4.0.5/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/extras/qd/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfsexports/SOFTWARE/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin
export vasp_std=/ibiscostorage/VirtualMatterLab/vasp_6.3.0_bin/vasp_std
mpirun -np 4 $vasp_std
