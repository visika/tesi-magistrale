import numpy as np
from ase.io import read
from mace.calculators.mace import MACECalculator

atoms=read('POSCAR','0') # define Atoms object by reading a geometry input file

path_to_model='.'
model_path=str(path_to_model)+'/MACE_ice13.model' # define complete path to the model

mace_ice=MACECalculator(str(model_path), default_dtype='float64',device='cuda') # define the calculator with the path to the model

atoms.calc=mace_ice # attach the calculator to the Atoms object
atoms.get_potential_energy() # compute the total energy
