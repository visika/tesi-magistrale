import numpy as np
import ase
from ase.io import read , write
import os


trj=read('data_for_train.extxyz',':')

for j in range(len(trj)):
	geo=trj[j].copy()
	os.mkdir(str(j+1))
	write(images=geo,filename=str(j+1)+'/POSCAR',format='vasp')
