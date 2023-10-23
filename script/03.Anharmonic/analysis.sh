#!/usr/bin/env sh

if [ ! -d analysis ]; then
    mkdir analysis
fi
# DONE Ricava U_0 da OUTCAR
U0=$(grep xlambda 01/OUTCAR | head -n 1 | awk '{print $6}')
echo "$U0"
# DONE Raccogli valori di OSZICAR
if [ -f analysis/df ]; then
    cat >analysis/df <<EOF
# xlambda <dF-U0>
EOF
fi
for folder in 01 02 03 04; do
    xlambda=$(grep xlambda $folder/OUTCAR | head -n 1 | awk '{print $3}')
    grep T= $folder/OSZICAR >"analysis/t.$xlambda"
    awk 'END{print NR}' "analysis/t.$xlambda"
    awk -v u="$U0" -v x="$xlambda" '{a+=$7+u}END{print x, a/NR/64}' "analysis/t.$xlambda" >>"analysis/df"
done
# DONE Calcola la media dell'energia di Helmholtz
# DONE Raccogli le medie per ogni lambda in un file
# TODO Integra per ottenere il volume corretto
