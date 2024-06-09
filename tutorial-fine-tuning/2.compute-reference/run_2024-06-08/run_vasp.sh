#!/usr/bin/env bash
poetry run python split_poscar.py
echo "Python finished, folders created"

number_of_structures=$(grep -c Lattice data_for_train.extxyz)

for ((j = 1; j <= number_of_structures; j++)); do
    cp INCAR KPOINTS POTCAR sbatch_template.sh $j

    cd $j || exit
    echo "Sbatching $j"
    sbatch sbatch_template.sh
    cd ..
done
