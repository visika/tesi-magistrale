#!/usr/bin/env python3
from ase.constraints import FixBondLength
from ase.io import read, write
import numpy as np
from ase.optimize import BFGS
from mace.calculators import mace_mp
import pandas as pd

atoms = read("init.xyz")

calc = mace_mp(
    model="/ibiscostorage/VirtualMatterLab/MACE-ICE13/MACE-ICE13-1.model",
    default_dtype="float64",
    device="cuda",
)

atoms_array = []
distance_oo = atoms.get_distance(0, 3)
vector_oo = atoms.get_distance(0, 3, mic=True, vector=True)
direction_oo = vector_oo / distance_oo

distances = np.linspace(2.0, 6.0, 15)

displacements = distances - distance_oo
potentialenergies = []

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
        logfile=f"optimization_dist={dist}.log",
        trajectory=f"optimization_dist={dist}.traj",
    )
    opt.run(steps=100)
    write(f"final_dist={dist}.xyz", atoms_distanced)
    potentialenergy = atoms_distanced.get_potential_energy()
    potentialenergies.append(potentialenergy)
    atoms_array.append(atoms_distanced)

dict = {"distance": distances, "energy": potentialenergies}
df = pd.DataFrame(data=dict)
df.to_csv("energies.csv", index=False)
