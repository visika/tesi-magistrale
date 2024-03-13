import numpy as np
from ase.io import read, write
import sys


file=sys.argv[1]

trj=read(str(file),index=':')

total_number_of_structures=len(trj)
number_of_selected_structures=int(total_number_of_structures*0.80)


train=[]
test=[]

for i in range(number_of_selected_structures):
    train.append(trj[i])

for i in range(number_of_selected_structures,total_number_of_structures):
    test.append(trj[i])




path='./'

for i in range(len(train)):
    write(images=train[i],filename=str(path)+'/training_set.extxyz',format='extxyz',append=True)

for i in range(len(test)):
    write(images=test[i],filename=str(path)+'/test_set.extxyz',format='extxyz',append=True)
