#!/bin/bash
# Get inside a poetry shell before running this

# Definisci il file di output
output_file="risultati.txt"

# Se il file di output esiste già, lo rimuove
# if [ -f "$output_file" ]; then
#     rm "$output_file"
# fi

# Loop attraverso i valori della variabile di input
for run in {1..3}; do
    echo "Run $run"
    # Esegue lo script Python con la variabile di input
    (./usr/bin/time --format="user+kernel %U+%S\nmax_memory %M" python script.py) > output.tmp 2>&1

    # Estrae le ultime tre linee dell'output
    ultime_tre_linee=$(cat output.tmp | tail -n 3)

    # Estrae i numeri dalle ultime tre linee
    num1=$(echo "$ultime_tre_linee" | sed -n '1p' | awk '{print $NF}')
    num2=$(echo "$ultime_tre_linee" | sed -n '2p' | awk '{print $NF}')
    num3=$(echo "$ultime_tre_linee" | sed -n '3p' | awk '{print $NF}')

    # Scrive i numeri nel file di output
    echo "$num1 $num2 $num3" >> "$output_file"
done
