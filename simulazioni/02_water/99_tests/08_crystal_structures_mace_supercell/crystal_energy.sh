#!/usr/bin/env sh
# Copia gli script in ciascuna cartella ed esegui separatamente i job di SLURM.
for poly in Ih II III IV VI VII VIII IX XI XIII XIV XV XVII
do
    cp crystal_energy.py $poly
    cd $poly
    sbatch ../slurm.sh
    cd ..
done
