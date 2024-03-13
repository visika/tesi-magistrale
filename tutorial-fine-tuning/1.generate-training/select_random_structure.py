import numpy as np
import ase 
from ase.io import read, write 
from ase.io.trajectory import Trajectory

import random

#trj=read('md.xyz', index=':')
trj=Trajectory('md.traj','r')

rnd = [random.randint(250, len(trj)-1) for _ in range(5)] # change 5 to how many structures you want
print(rnd)

rnd_geo=[]
for j in rnd:
    rnd_geo.append(trj[j])

write(images=rnd_geo, format='extxyz', filename='data_for_train.extxyz')
