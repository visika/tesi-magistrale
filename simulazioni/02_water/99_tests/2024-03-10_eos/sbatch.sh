#!/usr/bin/env sh
#SBATCH --job-name=eos
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

conda activate tesi_magistrale

echo Using "$(python -V)"
echo Python path: "$(which python)"

python3 ./eos.py --model medium --d3 --system structs/ICE13/Ih/POSCAR --device cuda

date
echo "FINISHED!!!"
rm $ALERT