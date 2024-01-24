#!/usr/bin/env python3
import os
from ase.io import read, write
from ase.optimize import BFGS
from mace.calculators import mace_mp
from ase.constraints import UnitCellFilter
from ase import units

atoms = read("POSCAR")
calc = mace_mp(model="medium", dispersion=False, default_dtype="float64", device="cuda")
atoms.calc = calc

path = "relax-full/"
os.makedirs(path, exist_ok=True)

for pressure in [0.0, 1.0]:
    geo = atoms.copy()
    geo.calc = calc
    ucf = UnitCellFilter(geo, scalar_pressure=pressure * units.GPa)
    opt = BFGS(
        ucf,
        logfile=path + f"optimization_{pressure}GPa.log",
        trajectory=path + f"optimization_{pressure}GPa.traj",
    )
    opt.run(fmax=0.001, steps=1000)
    write(images=geo, filename=path + f"final_{pressure}GPa.pdb")
