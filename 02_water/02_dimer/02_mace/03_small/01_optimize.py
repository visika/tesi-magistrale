#!/usr/bin/env python3
import os
from ase.io import write, read
from mace.calculators import mace_mp
from ase.optimize import BFGS

# Il file init.xyz contiene la geometria del dimero vicina a quella che
# desideriamo.
atoms = read("init.xyz")
atoms.center(vacuum=3.0)
calc = mace_mp(model="small", dispersion=False, default_dtype="float64", device="cuda")
atoms.calc = calc

# Percorso dove conservare le informazioni sul rilassamento.
path = "01_relax-positions/"
os.makedirs(path, exist_ok=True)

opt = BFGS(
    atoms,
    logfile=path + "optimization.log",
    trajectory=path + "optimization.traj",
)
opt.run(fmax=0.0001, steps=1000)
write(path + "final.pdb", atoms)