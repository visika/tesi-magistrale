#!/usr/bin/env python3
import os
from ase.io import read, write
from ase.optimize import BFGS
from mace.calculators import mace_mp

atoms = read("POSCAR")
calc = mace_mp(model="medium", dispersion=False, default_dtype="float64", device="cuda")
atoms.calc = calc

path = "relax-positions/"
os.makedirs(path, exist_ok=True)
opt = BFGS(
    atoms,
    logfile=path + "optimization.log",
    trajectory=path + "optimization.traj",
)
opt.run(fmax=0.001, steps=1000)
write(images=atoms, filename=path + "final.pdb")
