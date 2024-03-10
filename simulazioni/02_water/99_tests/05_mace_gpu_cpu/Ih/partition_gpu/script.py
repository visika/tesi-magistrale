#!/usr/bin/env python3
import os
from ase import units
from ase.constraints import UnitCellFilter
from ase.io import read, write
from ase.optimize import BFGS
from mace.calculators import mace_mp

atoms = read("POSCAR")
devices = ["cuda", "cpu"]

for device in devices:
    calc = mace_mp(
        model="medium", dispersion=False, default_dtype="float64", device=device
    )

    directory = f"relax-{device}"
    os.makedirs(directory, exist_ok=True)

    pressure = 0.0
    geo = atoms.copy()
    geo.calc = calc
    ucf = UnitCellFilter(geo, scalar_pressure=0.0 * units.GPa)
    opt = BFGS(
        ucf,
        logfile=f"{directory}/optimization.log",
        trajectory=f"{directory}/optimization.traj",
    )
    opt.run(fmax=0.0001, steps=1000)
    write(images=geo, filename=f"{directory}/final.pdb")

    nmol = geo.get_global_number_of_atoms() / 3.0
    potentialenergy = geo.get_potential_energy()
    e_crys = potentialenergy / nmol

    f = open(f"{directory}/e_crys.txt", "w")
    f.write(str(e_crys))
    f.close()
