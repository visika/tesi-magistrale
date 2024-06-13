python split_poscar.py

number_of_structures=`grep Lattice data_for_train.extxyz | wc -l`

for ((j=1; j<=$number_of_structures; j++))
do

cd $j/

cp	../INCAR .
cp 	../KPOINTS .
cp	../POTCAR .

cd ../
done
