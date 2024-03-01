#!/usr/bin/env python3
import os
from ase.io import read, write
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from mace.calculators import mace_mp

models = ["small", "medium", "large"]

for model in models:
    os.makedirs(model, exist_ok=True)
    os.chdir(model)
    atoms = read("../init.xyz")
    atoms.center(vacuum=50.0)
    calc = mace_mp(
        model=model, dispersion=False, default_dtype="float64", device="cuda"
    )
    atoms.calc = calc

    opt = BFGS(atoms, logfile="optimization.log", trajectory="optimization.traj")
    opt.run(fmax=0.0001, steps=1000)
    write("final.xyz", atoms)

    for d in [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]:
        vib = Vibrations(atoms, delta=d, name=f"vib_delta={d}")
        vib.run()
        vib.summary(log=f"dimer_delta={d}_summary.txt")
        vib.clean()

    os.chdir("..")
