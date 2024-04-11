#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
from ase.io import read, write
from mace.calculators import mace_mp
from ase.optimize import BFGS
from ase.constraints import FixBondLength

# Leggi la geometria iniziale del dimero
atoms = read("init.xyz")

calc = mace_mp(model="medium", dispersion=True, default_dtype="float64", device="cuda")

# Costruisco un array nel quale conservare tutte le configurazioni di atomi
# posti a distanze differenti tra loro
atoms_array = []
distance_oo = atoms.get_distance(0, 3)
vector_oo = atoms.get_distance(0, 3, mic=True, vector=True)
direction_oo = vector_oo / distance_oo
# Definisco l'intervallo di distanze tra i due atomi di ossigeno che voglio
# campionare. Qualcosa come un punto ogni 0.1 Å può andare bene.
distances = np.linspace(2.0, 6.0, num=40)
displacements = distances - distance_oo
potentialenergies = []

os.makedirs("geometries", exist_ok=True)
os.makedirs("optimization", exist_ok=True)
for disp, dist in zip(displacements, distances):
    atoms_distanced = atoms.copy()
    atoms_distanced[3].position += direction_oo * disp
    atoms_distanced[4].position += direction_oo * disp
    atoms_distanced[5].position += direction_oo * disp
    constraint = FixBondLength(0, 3)
    atoms_distanced.set_constraint(constraint)
    atoms_distanced.calc = calc
    opt = BFGS(
        atoms_distanced,
        logfile=f"optimization/dist={dist}.log",
        trajectory=f"optimization/dist={dist}.traj",
    )
    opt.run(steps=100)
    write(f"geometries/final_dist={dist}.xyz", atoms_distanced)
    potentialenergy = atoms_distanced.get_potential_energy()
    potentialenergies.append(potentialenergy)
    atoms_array.append(atoms_distanced)

dict = {"distance": distances, "energy": potentialenergies}
df = pd.DataFrame(data=dict)
df.to_csv("energies.csv", index=False)
