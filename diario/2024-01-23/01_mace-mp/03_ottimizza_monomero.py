#!/usr/bin/env python3
from mace.calculators import mace_mp
from ase import build
from ase.optimize import BFGS
import os
from ase.io import write

atoms = build.molecule("H2O")
calc = mace_mp(model="medium", dispersion=False, default_dtype="float64", device="cuda")
atoms.calc = calc
print("Potential energy:", atoms.get_potential_energy())

path = "relax-monomer/"
os.makedirs(path, exist_ok=True)

opt = BFGS(
    atoms, trajectory=path + "optimization.traj", logfile=path + "optimization.log"
)
opt.run(fmax=0.001)

write(images=atoms, filename=path + "final.pdb")
