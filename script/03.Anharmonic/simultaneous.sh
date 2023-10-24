#!/usr/bin/env bash
# DONE Adatta lo script al sistema ibisco
# Procedura dalla lezione del 2023-05-12
if [ ! -f POTCAR ]; then
   cp ../00.PerfectCrystal/POTCAR .
fi
# Vogliamo una super-cella con 64 atomi, cioè 4x4x4
# Scegliamo un volume di 16.6 * 64 = 1062.4
s=4
k=1
v=16.6
vsuper=1062.4
t=600

cat >POSCAR <<EOF
super cell
        -$vsuper
   2.000000000000000   0.000000000000000   2.000000000000000
   2.000000000000000   2.000000000000000   0.000000000000000
   0.000000000000000   2.000000000000000   2.000000000000000
   64
 Direct
   0.000000000000000   0.000000000000000   0.000000000000000
   0.250000000000000   0.000000000000000   0.000000000000000
   0.500000000000000   0.000000000000000   0.000000000000000
   0.750000000000000   0.000000000000000   0.000000000000000
   0.000000000000000   0.250000000000000   0.000000000000000
   0.250000000000000   0.250000000000000   0.000000000000000
   0.500000000000000   0.250000000000000   0.000000000000000
   0.750000000000000   0.250000000000000   0.000000000000000
   0.000000000000000   0.500000000000000   0.000000000000000
   0.250000000000000   0.500000000000000   0.000000000000000
   0.500000000000000   0.500000000000000   0.000000000000000
   0.750000000000000   0.500000000000000   0.000000000000000
   0.000000000000000   0.750000000000000   0.000000000000000
   0.250000000000000   0.750000000000000   0.000000000000000
   0.500000000000000   0.750000000000000   0.000000000000000
   0.750000000000000   0.750000000000000   0.000000000000000
   0.000000000000000   0.000000000000000   0.250000000000000
   0.250000000000000   0.000000000000000   0.250000000000000
   0.500000000000000   0.000000000000000   0.250000000000000
   0.750000000000000   0.000000000000000   0.250000000000000
   0.000000000000000   0.250000000000000   0.250000000000000
   0.250000000000000   0.250000000000000   0.250000000000000
   0.500000000000000   0.250000000000000   0.250000000000000
   0.750000000000000   0.250000000000000   0.250000000000000
   0.000000000000000   0.500000000000000   0.250000000000000
   0.250000000000000   0.500000000000000   0.250000000000000
   0.500000000000000   0.500000000000000   0.250000000000000
   0.750000000000000   0.500000000000000   0.250000000000000
   0.000000000000000   0.750000000000000   0.250000000000000
   0.250000000000000   0.750000000000000   0.250000000000000
   0.500000000000000   0.750000000000000   0.250000000000000
   0.750000000000000   0.750000000000000   0.250000000000000
   0.000000000000000   0.000000000000000   0.500000000000000
   0.250000000000000   0.000000000000000   0.500000000000000
   0.500000000000000   0.000000000000000   0.500000000000000
   0.750000000000000   0.000000000000000   0.500000000000000
   0.000000000000000   0.250000000000000   0.500000000000000
   0.250000000000000   0.250000000000000   0.500000000000000
   0.500000000000000   0.250000000000000   0.500000000000000
   0.750000000000000   0.250000000000000   0.500000000000000
   0.000000000000000   0.500000000000000   0.500000000000000
   0.250000000000000   0.500000000000000   0.500000000000000
   0.500000000000000   0.500000000000000   0.500000000000000
   0.750000000000000   0.500000000000000   0.500000000000000
   0.000000000000000   0.750000000000000   0.500000000000000
   0.250000000000000   0.750000000000000   0.500000000000000
   0.500000000000000   0.750000000000000   0.500000000000000
   0.750000000000000   0.750000000000000   0.500000000000000
   0.000000000000000   0.000000000000000   0.750000000000000
   0.250000000000000   0.000000000000000   0.750000000000000
   0.500000000000000   0.000000000000000   0.750000000000000
   0.750000000000000   0.000000000000000   0.750000000000000
   0.000000000000000   0.250000000000000   0.750000000000000
   0.250000000000000   0.250000000000000   0.750000000000000
   0.500000000000000   0.250000000000000   0.750000000000000
   0.750000000000000   0.250000000000000   0.750000000000000
   0.000000000000000   0.500000000000000   0.750000000000000
   0.250000000000000   0.500000000000000   0.750000000000000
   0.500000000000000   0.500000000000000   0.750000000000000
   0.750000000000000   0.500000000000000   0.750000000000000
   0.000000000000000   0.750000000000000   0.750000000000000
   0.250000000000000   0.750000000000000   0.750000000000000
   0.500000000000000   0.750000000000000   0.750000000000000
   0.750000000000000   0.750000000000000   0.750000000000000
EOF

cat >SPOSCAR <<EOF
super cell
        4.0000000000
   2.000000000000000   0.000000000000000   2.000000000000000
   2.000000000000000   2.000000000000000   0.000000000000000
   0.000000000000000   2.000000000000000   2.000000000000000
   64
 Direct
   0.000000000000000   0.000000000000000   0.000000000000000
   0.250000000000000   0.000000000000000   0.000000000000000
   0.500000000000000   0.000000000000000   0.000000000000000
   0.750000000000000   0.000000000000000   0.000000000000000
   0.000000000000000   0.250000000000000   0.000000000000000
   0.250000000000000   0.250000000000000   0.000000000000000
   0.500000000000000   0.250000000000000   0.000000000000000
   0.750000000000000   0.250000000000000   0.000000000000000
   0.000000000000000   0.500000000000000   0.000000000000000
   0.250000000000000   0.500000000000000   0.000000000000000
   0.500000000000000   0.500000000000000   0.000000000000000
   0.750000000000000   0.500000000000000   0.000000000000000
   0.000000000000000   0.750000000000000   0.000000000000000
   0.250000000000000   0.750000000000000   0.000000000000000
   0.500000000000000   0.750000000000000   0.000000000000000
   0.750000000000000   0.750000000000000   0.000000000000000
   0.000000000000000   0.000000000000000   0.250000000000000
   0.250000000000000   0.000000000000000   0.250000000000000
   0.500000000000000   0.000000000000000   0.250000000000000
   0.750000000000000   0.000000000000000   0.250000000000000
   0.000000000000000   0.250000000000000   0.250000000000000
   0.250000000000000   0.250000000000000   0.250000000000000
   0.500000000000000   0.250000000000000   0.250000000000000
   0.750000000000000   0.250000000000000   0.250000000000000
   0.000000000000000   0.500000000000000   0.250000000000000
   0.250000000000000   0.500000000000000   0.250000000000000
   0.500000000000000   0.500000000000000   0.250000000000000
   0.750000000000000   0.500000000000000   0.250000000000000
   0.000000000000000   0.750000000000000   0.250000000000000
   0.250000000000000   0.750000000000000   0.250000000000000
   0.500000000000000   0.750000000000000   0.250000000000000
   0.750000000000000   0.750000000000000   0.250000000000000
   0.000000000000000   0.000000000000000   0.500000000000000
   0.250000000000000   0.000000000000000   0.500000000000000
   0.500000000000000   0.000000000000000   0.500000000000000
   0.750000000000000   0.000000000000000   0.500000000000000
   0.000000000000000   0.250000000000000   0.500000000000000
   0.250000000000000   0.250000000000000   0.500000000000000
   0.500000000000000   0.250000000000000   0.500000000000000
   0.750000000000000   0.250000000000000   0.500000000000000
   0.000000000000000   0.500000000000000   0.500000000000000
   0.250000000000000   0.500000000000000   0.500000000000000
   0.500000000000000   0.500000000000000   0.500000000000000
   0.750000000000000   0.500000000000000   0.500000000000000
   0.000000000000000   0.750000000000000   0.500000000000000
   0.250000000000000   0.750000000000000   0.500000000000000
   0.500000000000000   0.750000000000000   0.500000000000000
   0.750000000000000   0.750000000000000   0.500000000000000
   0.000000000000000   0.000000000000000   0.750000000000000
   0.250000000000000   0.000000000000000   0.750000000000000
   0.500000000000000   0.000000000000000   0.750000000000000
   0.750000000000000   0.000000000000000   0.750000000000000
   0.000000000000000   0.250000000000000   0.750000000000000
   0.250000000000000   0.250000000000000   0.750000000000000
   0.500000000000000   0.250000000000000   0.750000000000000
   0.750000000000000   0.250000000000000   0.750000000000000
   0.000000000000000   0.500000000000000   0.750000000000000
   0.250000000000000   0.500000000000000   0.750000000000000
   0.500000000000000   0.500000000000000   0.750000000000000
   0.750000000000000   0.500000000000000   0.750000000000000
   0.000000000000000   0.750000000000000   0.750000000000000
   0.250000000000000   0.750000000000000   0.750000000000000
   0.500000000000000   0.750000000000000   0.750000000000000
   0.750000000000000   0.750000000000000   0.750000000000000
EOF

for folder in 00 01 02 03 04 05; do
   mkdir $folder
   cp POSCAR $folder
done

cat >01/INCAR <<EOF
XLAMBDA=0.0
EOF

cat >02/INCAR <<EOF
XLAMBDA=0.333333333
EOF

cat >03/INCAR <<EOF
XLAMBDA=0.666666666
EOF

cat >04/INCAR <<EOF
XLAMBDA=1.0
EOF

# Allo stato dei fatti, abbiamo una struttura con tutti gli atomi
# nelle proprie posizioni di equilibrio. Avviando VASP con questa
# configurazione, otterremmo un valore dell'energia corrispondente a
# U_0. Dobbiamo configurare INCAR per abilitare la dinamica
# molecolare. Il parametro per attivare la dinamica molecolare è
# IBRION=0.

# TODO NSW dovrebbe essere 2000 per avere buoni risultati, ma ora sto
# facendo un test, quindi lo imposto a 200.

# POTIM=3.0 è il passo temporale in femto secondi. TEBEG=600 è la
# temperatura iniziale della simulazione.

# Usiamo EDIFF=1d-4 perché la convergenza sia meno stringente e più
# rapida. IALGO=48 impiega un algoritmo di diagonalizzazione più
# rapido ma anche più instabile. MAXMIX=40 utilizza più densità dai
# passi precedenti, per velocizzare la convergenza.

# TIPO=inteharm indica che vogliamo effettuare l'integrazione
# termodinamica con una combinazione di potenziale armonico e
# potenziale completo.

cat >INCAR.molecular-dynamics <<EOF
EDIFF=1d-4
SIGMA=0.544
LWAVE=.FALSE.
LCHARG=.FALSE.
IMAGES=4
SPRING=-1000

#ENCUT=400

NPAR=1

IALGO=48
IBRION=0
NSW=500
#ISYM=0
POTIM=3.0
TEBEG=$t
LANDERSON=.T.
NANDERSON=50
#POMASS=26.98
#NBLOCK=5
#NWRITE=0
MAXMIX=40

TIPO=inteharm
#TIPO=harmonic
XMIXHARM=1.0
XMIXREF=0.0
R1=8.05
NCELL=64

XLAMBDA=0.0
#LADIABATIC=.T.
EOF

# Dato che ci interessano solo le differenze rispetto al valore di
# riferimento dell'energia, effettuiamo i calcoli solo in un punto.

cat >KPOINTS <<EOF
auto
 0
Gamma
1 1 1
0 0 0
EOF

# We need the FORCES file to generate the HARMONIC file
# TODO Crea flusso di creazione del file con runphon

cat >INCAR.harmonic <<EOF
EDIFF=1d-6
SIGMA=0.544
LWAVE=.FALSE.
LCHARG=.FALSE.
#ENCUT=400
NPAR=1
EOF

FILE_FORCES=FORCES.s$s.k$k.v$v
if [ -f FORCES ]; then
  rm FORCES
fi
if [ -f "$SLURM_SUBMIT_DIR/$FILE_FORCES" ]; then
  echo "File $FILE_FORCES found"
  ln -Ts "$SLURM_SUBMIT_DIR/$FILE_FORCES" FORCES
else
   echo "Running runphon to generate $FILE_FORCES"
   cp INCAR.harmonic INCAR
   $RUNPHON
   cp FORCES "$SLURM_SUBMIT_DIR/$FILE_FORCES"
fi
unset FILE_FORCES

# Remember to use the correct HARMONIC file force constant matrix and set the
# correct radius TODO Creare il file HARMONIC adatto al volume, se non esiste
# già Per far scrivere a PHON nel file HARMONIC la matrice delle costanti di
# forza, impostare nel file INPHON le seguenti variabili:
# LSUPER = .F., LFORCEOUT = .T.

if [ ! -f HARMONIC.v$v ]; then
   echo "    HARMONIC.v$v not found; generating it with PHON"
   cat >INPHON <<EOF
# number of ions types and masses
NTYPES = 1
MASS =26.98
LCENTRAL=.FALSE.
#LSYMM=.F.

# generate superlattice
LSUPER = .FALSE.
NDIM = 4 4 4
DXSTART = 1 -1 1
LFORCEOUT=.TRUE.

# free energy calculation
LFREE = .FALSE.
TEMPERATURE = 29
#PTEMP = 10 81

# q points section
LRECIP = .FALSE.
ND = 3
NPOINTS = 51
QI = 0.0        0.0        0.0      1.0        1.0        0.0      0.0        0.0        0.0

QF = 1.0        0.0        0.0      0.0        0.0        0.0      0.5        0.5        0.5

# density of states
LGAMMA = .FALSE.
QA = 8; QB = 8; QC = 8
DOSIN = 0
DOSEND = 10
DOSSTEP = 0.05
DOSSMEAR = 0.05

# verbosity
IPRINT = 0
EOF
   $PHON
fi

cp HARMONIC.v$v HARMONIC

# Run VASP
cp INCAR.molecular-dynamics INCAR
$MPI -np 4 "$VASPGAMMA"
