#!/usr/bin/env python3
import os
from ase.io import read
from mace.calculators import mace_mp
from ase.vibrations import Vibrations

path = "01_relax-positions/"
atoms = read(path + "final.pdb")
atoms.center(vacuum=3.0)
calc = mace_mp(model="medium", dispersion=False, default_dtype="float64", device="cuda")
atoms.calc = calc

path = "02_vib/"
os.makedirs(path, exist_ok=True)
for d in [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]:
    vib = Vibrations(atoms, delta=d, name=f"vib_delta={d}")
    vib.run()
    vib.summary(log=path + f"dimer_delta={d}_summary.txt")
