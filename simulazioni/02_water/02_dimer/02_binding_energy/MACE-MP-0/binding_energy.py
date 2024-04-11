#!/usr/bin/env python3
import numpy as np
import pandas as pd
from ase.io import read, write
from mace.calculators import mace_mp
from ase.optimize import BFGS
from ase.constraints import FixBondLength

# Leggi la geometria iniziale del dimero
atoms = read("init.xyz")

calc = mace_mp(model="medium", dispersion=True, default_dtype="float64", device="cuda")

atoms_array = []
distance_oo = atoms.get_distance(0, 3)
vector_oo = atoms.get_distance(0, 3, mic=True, vector=True)
direction_oo = vector_oo / distance_oo
distances = np.linspace(2.0, 2.7, 8)
displacements = distances - distance_oo
potentialenergies = []

for d, dist in zip(displacements, distances):
    atoms_distanced = atoms.copy()
    atoms_distanced[3].position += direction_oo * d
    atoms_distanced[4].position += direction_oo * d
    atoms_distanced[5].position += direction_oo * d
    constraint = FixBondLength(0, 3)
    atoms_distanced.set_constraint(constraint)
    atoms_distanced.calc = calc
    opt = BFGS(
        atoms_distanced,
        logfile=f"optimization_d={d}.log",
        trajectory=f"optimization_d={d}.traj",
    )
    opt.run(steps=100)
    write(f"final_d={dist}.xyz", atoms_distanced)
    potentialenergy = atoms_distanced.get_potential_energy()
    potentialenergies.append(potentialenergy)
    atoms_array.append(atoms_distanced)

dict = {"distance": distances, "energy": potentialenergies}
df = pd.DataFrame(data=dict)
df.to_csv("energies_2.csv", index=False)
