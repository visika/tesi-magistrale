#!/usr/bin/env sh
# Questo programma genera una super-cella con dimensione specificata dalla variabile -s passata in input.
# Il file 02_harmonic.sh calcola poi le forze e produce in output le frequenze fononiche.

while getopts s:k:d: flag
do
     case "${flag}" in
          s) s=${OPTARG};;
          k) k=${OPTARG};;
          d) d=${OPTARG};;
     esac
done

if [ -z $s ]; then
    echo "Variable s not defined, stopping now"
    exit 1
fi

if [ -z $k ]; then
    echo "Variable k not defined, stopping now"
    exit 1
fi

echo "Supercell size: $s";
echo "K-points sampling: $k";
echo "Custom displacement: $d";

# cp /data/studenti/mollo/aluminium/POTCAR.bak POTCAR

# Il file POSCAR viene sovrascritto con i contenuti di SPOSCAR.
# È bene reinizializzarlo ogni volta.
cat > POSCAR<<EOF
Ice Ih
1.0
  7.6780930000 0.0000000000 0.0000000000 
  3.8390460000 6.6494230000 0.0000000000 
  0.0000000000 0.0000000000 7.2345670000 
H O
24 12
Direct
 0.000007000  0.334718000  0.199799000
 0.665252000  0.000010000  0.199780000
 0.334714000  0.665261000  0.199791000
 0.334760000  0.999976000  0.699786000
 0.999980000  0.665239000  0.699800000
 0.665240000  0.334755000  0.699799000
 0.544461000  0.000006000  0.019584000
-0.000011000  0.455520000  0.019605000
 0.455505000  0.544480000  0.019594000
 0.455543000  0.999980000  0.519584000
 0.999998000  0.544446000  0.519600000
 0.544461000  0.455524000  0.519592000
 0.332240000  0.879481000  0.984486000
 0.211731000  1.120507000  0.984481000
 0.879462000  0.788265000  0.984489000
 0.788248000  0.332255000  0.984498000
 0.667731000  0.211748000  0.984486000
 0.120526000  0.211703000  0.484497000
 0.667762000  1.120506000  0.484488000
 0.788272000  0.879477000  0.484480000
 0.211722000  0.667743000  0.484495000
 0.332240000  0.788252000  0.484486000
-0.120503000  0.332222000  0.484499000
 1.120486000  0.667750000  0.984491000
 0.000009000  0.331330000  0.061616000
 0.668637000  0.000012000  0.061595000
 0.331326000  0.668660000  0.061608000
 0.331370000  0.999972000  0.561601000
 0.668637000  0.331348000  0.561616000
 0.335676000  0.999992000  0.936746000
 0.664294000  0.335703000  0.936768000
 0.000012000  0.335655000  0.436771000
 0.664333000  0.999995000  0.436740000
 0.335670000  0.664304000  0.436760000
 0.999979000  0.668634000  0.561617000
 0.999974000  0.664306000  0.942299010
EOF

cat > INCAR <<EOF
NPAR=4
ISMEAR=0
SIGMA=0.01
IVDW=11
LWAVE=.FALSE.
LCHARG=.FALSE.
LREAL=Auto
EOF

# Generazione della supercella
cat > INPHON <<EOF
# Generated by harmonic.sh

# number of ions types and masses
NTYPES = 2

# generate superlattice
LSUPER = .TRUE.
NDIM = $s $s $s
EOF

cat > KPOINTS <<EOF
auto
 0
Gamma
$k $k $k
0 0 0
EOF

# Run phon to get the file SPOSCAR
rm POSCAR.phon
echo "Running PHON to generate the supercell"
./phon
# SPOSCAR and DISP are changed

# Copy SPOSCAR to POSCAR
cp SPOSCAR POSCAR

echo "Now put the contents of DISP inside runphon, then run harmonic_2.sh"

# if [ -z $d ]; then
#     # We need to change the suggested displacement in runphon
#     export d=$(cat DISP | rev | cut -c 2- | rev)
#     sed "24s/.*/$d/" /data/studenti/mollo/aluminium/runphon.bak > runphon
# else
#     export x="   1  $d -$d  $d "
#     sed "24s/.*/\x22$x\x22/" /data/studenti/mollo/aluminium/runphon.bak > runphon
# fi

exit 0
