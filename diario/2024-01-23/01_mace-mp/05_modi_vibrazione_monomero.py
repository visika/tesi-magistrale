#!/usr/bin/env python3
from ase.vibrations.vibrations import Vibrations
from mace.calculators import mace_mp
from ase import Atoms
from ase.optimize import BFGS
import os
from ase.io import read, write

# atoms = build.molecule("H2O")
# atoms = Atoms("H2O", positions=[(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)])
atoms = read("relax-monomer-hard/final.pdb")
# atoms.center(vacuum=3.0)
calc = mace_mp(model="medium", dispersion=False, default_dtype="float64", device="cuda")
atoms.calc = calc
# print("Potential energy:", atoms.get_potential_energy())

# path = "relax-monomer-hard/"
# os.makedirs(path, exist_ok=True)

# opt = BFGS(
# atoms, trajectory=path + "optimization.traj", logfile=path + "optimization.log"
# )
# opt.run(fmax=0.001)

# write(images=atoms, filename=path + "final.pdb")

vib = Vibrations(atoms)
vib.run()
vib.summary(log="H2O_MACE_summary.txt")
vib.write_mode()
