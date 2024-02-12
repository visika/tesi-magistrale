#!/usr/bin/env bash
#SBATCH --job-name=H2O
#SBATCH --nodes=1             # Number of nodes
#SBATCH --tasks-per-node=1    # Number of MPI ranks per node
#SBATCH --cpus-per-task=4
#SBATCH --partition=parallel

ulimit -s unlimited

# export OMP_NUM_THREADS=2
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

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/lustre/home/mmollo/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/lustre/home/mmollo/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/lustre/home/mmollo/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/lustre/home/mmollo/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate mace_env

echo Using "$(python -V)"
python 01_optimize.py

echo "FINISHED!!!"
rm $ALERT
