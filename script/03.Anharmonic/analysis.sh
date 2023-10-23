#!/usr/bin/env sh

if [ ! -d analysis ]; then
    mkdir analysis
fi
# DONE Ricava U_0 da OUTCAR
U0=$(grep xlambda 01/OUTCAR | head -n 1 | awk '{print $6}')
echo "$U0"
# TODO Raccogli valori di OSZICAR
for folder in 01 02 03 04; do
    xlambda=$(grep xlambda $folder/OUTCAR | head -n 1 | awk '{print $3}')
    grep T= $folder/OSZICAR >"analysis/t.$xlambda"
    awk 'END{print NR}' "analysis/t.$xlambda"
    awk -v u="$U0" '{a+=$7+u}END{print a/NR/64}' "analysis/t.$xlambda" >"analysis/l.$xlambda"
done
cat analysis/l* >analysis/df
# TODO Calcola la media dell'energia di Helmholtz
# TODO Raccogli le medie per ogni lambda in un file
# TODO Integra per ottenere il volume corretto
