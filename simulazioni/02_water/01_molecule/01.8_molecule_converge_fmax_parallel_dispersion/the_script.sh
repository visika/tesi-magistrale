#!/usr/bin/env bash
for model in small medium large; do
    for fmax in 1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8; do
        mkdir -p $model/$fmax
        cp sbatch.sh $model/$fmax
        cp converge_fmax.py $model/$fmax
        cd $model/$fmax || exit
        sed -i "s/mmm/$model/g" sbatch.sh
        sed -i "s/fff/$fmax/g" sbatch.sh
        sed -i "s/pfmax/$fmax/g" sbatch.sh
        sbatch sbatch.sh
        cd ../..
    done
done
