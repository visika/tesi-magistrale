#!/usr/bin/bash
#SBATCH --job-name=Mollo
#SBATCH --nodes=1              # Number of nodes
#SBATCH --tasks-per-node=1     # Number of MPI ranks per node
##SBATCH --exclusive
##SBATCH --exclude=ibiscohpc-wn04,ibiscohpc-wn06,ibiscohpc-wn07
#SBATCH --cpus-per-task=1
##SBATCH --mem=1400000
#SBATCH --partition=sequential

# Requires: POTCAR, fit, in

ulimit -s unlimited

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

source /nfsexports/SOFTWARE/intel/oneapi/setvars.sh
export MPI="mpirun"
echo "     using $(which $MPI)"

# Prendi il percorso di VASP
export VASP="/lustre/home/tccourse/vasp46-da/vasp"

echo "Il mio programma sta iniziando..."

# Passa alla cartella di esecuzione dello script
# Molto probabilmente è quella attuale
cd "$SLURM_SUBMIT_DIR" || exit

if [ -f vol.out ]; then
    rm vol.out
fi

cat >INCAR <<EOF
ENCUT=400
NPAR=4
EOF

for k in 1 2 4 6 8 10 12 14 16 18 20; do
    cat >KPOINTS <<EOF
auto
 0
Gamma
$k $k $k
0 0 0
EOF

    if [ -f u$k ]; then
        rm u$k
    fi
    for v in 15.0 15.2 15.4 15.6 15.8 16.0 16.2 16.4 16.6 16.8 17.0; do
        cat >POSCAR <<EOF
Aluminium crystal
        -$v
   0.500000000000000   0.000000000000000   0.500000000000000
   0.500000000000000   0.500000000000000   0.000000000000000
   0.000000000000000   0.500000000000000   0.500000000000000
    1
 Direct
   0.000000000000000   0.000000000000000   0.000000000000000
EOF

        $MPI $VASP >out

        grep "free  energy" OUTCAR | awk -v vol=$v '{print vol,$5}' >>u$k
    done
    cp u$k tmp
    ./fit <in | grep "v0 " | awk -v k=$k '{print k,$3}' >>vol.out
    rm tmp
done
