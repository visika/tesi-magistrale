import numpy as np
from ase.io import read
#from mace.calculators.mace import MACECalculator
from ase.md.langevin import Langevin # NVT molecular dynamics
from ase.md import MDLogger
import ase.units as units
from mace.calculators import mace_mp

#mace=MACECalculator('/rds/user/fd385/hpc-work/work/MolCrys/ICE/it_0/models/MACE_model_large_swa.model', default_dtype='float64',device='cuda')
mace=mace_mp()
atoms=read('POSCAR',0)
atoms.calc = mace


dyn = Langevin(
    atoms,
    timestep=0.5 * units.fs,
    temperature_K=100.0,  # temperature in K
    friction=0.01 / units.fs,
    trajectory='md.traj',
)

dyn.attach(MDLogger(dyn, atoms, 'md.log', header=True, stress=False, mode="a"), interval=4)
dyn.run(10000)
