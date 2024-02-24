#!/usr/bin/env sh
for poly in Ih II III IV VI VII VIII IX XI XIII XIV XV XVII
do
    cp crystal_energy.py $poly
    cd $poly
    sbatch ../slurm.sh
    cd ..
done
