import numpy as np
from ase.io import read
from mace.calculators.mace import MACECalculator

atoms=read('test.extxyz',':')

model_paths=[]
model_paths.append('./MACE_model_swa.model')

models=[]
models.append(MACECalculator(model_paths[0], default_dtype='float64',device='cuda'))

energy=[]
forces=[]
stress=[]
for i in range(len(atoms)):
    copy=atoms[i].copy()
    copy.calc=models[0]
    
    energy.append((atoms[i].get_potential_energy() - copy.get_potential_energy())/atoms[i].get_global_number_of_atoms())
    forces.append((atoms[i].get_forces() - copy.get_forces()))
    stress.append((atoms[i].get_stress() - copy.get_stress()))

energy=np.array(energy)*1000.
forces=np.concatenate(forces)*1000
stress=np.array(stress)*1000.

print('TEST SET ERROR')
print('RMSE Energy meV/atom: ',np.sqrt(np.mean(energy**2))) 
print('RMSE Forces meV/A: ',np.sqrt(np.mean(forces**2))) 
print('RMSE Stress meV/A3: ',np.sqrt(np.mean(stress**2))) 
