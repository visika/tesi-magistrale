#!/usr/bin/env sh
for model in small medium large; do
    mkdir $model
    cp sbatch.sh $model
    cp optimize.py $model
    cd $model
    sed -i "s/mmm/$model/g" sbatch.sh
    sbatch sbatch.sh
    cd ..
done
