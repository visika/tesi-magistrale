#!/usr/bin/env bash
#SBATCH --job-name=Mollo
#SBATCH --nodes=1               # Number of nodes
#SBATCH --tasks-per-node=1     # Number of MPI ranks per node
##SBATCH --exclusive
##SBATCH --exclude=ibiscohpc-wn04,ibiscohpc-wn06,ibiscohpc-wn07
#SBATCH --cpus-per-task=1
##SBATCH --mem=1400000
#SBATCH --partition=sequential

# Da eseguire nella cartella 03.Anharmonic

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

ALERT=000-ANALYSIS_IS_RUNNING
touch $ALERT

t=600

cat >t$t/anharmonic_correction <<EOF
# volume correction
EOF

cat >t$t/free.harmonic.tmp <<EOF
# volume helmholtz_energy
EOF

for v in 16.6 16.8 17.0; do
    echo "[I] Inizia loop di analisi per v=$v"
    cd "$SLURM_SUBMIT_DIR/t$t/v$v" || exit

    if [ ! -d analysis ]; then
        mkdir analysis
    fi

    # DONE Ricava U_0 da OUTCAR
    U0=$(grep xlambda 01/OUTCAR | head -n 1 | awk '{print $6}')
    echo "[I] L'energia di punto zero è $U0"

    # DONE Raccogli valori di OSZICAR
    # DONE Calcola la media dell'energia di Helmholtz
    # DONE Raccogli le medie per ogni lambda in un file

    cat >analysis/df <<EOF
# xlambda <DF+U0>
EOF

    for folder in 01 02 03 04; do
        xlambda=$(grep xlambda $folder/OUTCAR | head -n 1 | awk '{print $3}')
        grep T= $folder/OSZICAR >"analysis/t.$xlambda"
        echo "[I] Ci sono $(awk 'END{print NR}' "analysis/t.$xlambda") record di dinamica molecolare"
        awk -v u="$U0" -v x="$xlambda" '{a+=$7+u}END{print x, a/NR/64}' "analysis/t.$xlambda" >>"analysis/df"
        rm "analysis/t.$xlambda"
    done
    # TODO Integra per ottenere il volume corretto
    cd analysis || exit
    if [ -L integrate.py ] || [ -e integrate.py ]; then
        rm integrate.py
    fi
    ln -s "$SLURM_SUBMIT_DIR/integrate.py" .
    echo "[I] Calculating thermodynamic integral"
    python3 integrate.py >integral
    echo "[I] L'integrale per v=$v vale $(cat integral)"

    awk -v vol=v$v '{print vol, $1}' integral >>"$SLURM_SUBMIT_DIR/t$t/anharmonic_correction"
    # Prendi solo i valori di volume da confrontare dalle energie armoniche calcolate in precedenza
    grep $v "$SLURM_SUBMIT_DIR/t$t/free.s4.k2.t600" >>"$SLURM_SUBMIT_DIR/t$t/free.harmonic.tmp"

    cd "$SLURM_SUBMIT_DIR" || exit
done # v

# Alla fine di questo loop, ho un file integral con il valore dell'integrale
# termodinamico per ciascun volume. Il passo successivo è accorpare questi dati
# e confrontarli con l'approssimazione armonica.

cd "$SLURM_SUBMIT_DIR/t$t" || exit

# --lines=+2 ignora la prima riga, che è di intestazione per entrambi i file in ingresso
paste free.harmonic.tmp anharmonic_correction | tail --lines=+2 | awk '{print $1, $2+$4}' >helmholtz.plus
paste free.harmonic.tmp anharmonic_correction | tail --lines=+2 | awk '{print $1, $2-$4}' >helmholtz.minus

echo "[I] Nel file helmholtz trovi i valori dell'energia di Helmholtz per i volumi forniti in ingresso"
echo "[I]    FINISHED!!!"

cd "$SLURM_SUBMIT_DIR" || exit
rm $ALERT
