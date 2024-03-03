#!/usr/bin/env sh
#SBATCH --job-name=1e-6
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --propagate=NONE

ulimit -s unlimited

source /nfsexports/SOFTWARE/intel/oneapi/setvars.sh
export MPI="mpirun"
echo "     using $(which $MPI)"

cd "$SLURM_SUBMIT_DIR" || exit

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
python converge_fmax.py --model large --fmax 1e-6

echo "FINISHED!!!"
rm $ALERT
