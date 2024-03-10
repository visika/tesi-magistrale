#!/usr/bin/env sh
# Assicurarsi che esistano e che siano effettivamente directory gli elementi
# iterati con poly.
for poly in Ih II III IV VI VII VIII IX XI XIII XIV XV XVII
do
    cp crystal_energy.py $poly
    cd $poly
    sbatch ../slurm.sh
    cd ..
done
