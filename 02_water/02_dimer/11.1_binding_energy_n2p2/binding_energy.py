#!/usr/bin/env python3
from ase.constraints import FixBondLength
from ase.io import read, write
import numpy as np
from ase.optimize import BFGS
import pandas as pd
from ase.calculators.lammpslib import LAMMPSlib

# Inizializza LAMMPS
atom_types = {"O": 2, "H": 1}
cmd_nn = [
    "#!/bin/bash",
    'variable runnerDir       string "/ibiscostorage/VirtualMatterLab/train_008"',
    "variable runnerCutoff    equal  8.466835984",
    "pair_style hdnnp ${runnerCutoff} dir ${runnerDir} showew no showewsum 1 resetew yes maxew 200000 cflength 1.889726 cfenergy 0.036749",
    "pair_coeff * * H O",
]
calc = LAMMPSlib(
    atom_types=atom_types,
    lmpcmds=cmd_nn,
    tmp_dir="./",
    log_file="log.lammps",
    keep_alive=True,
)

atoms = read("init.xyz")

atoms_array = []
distance_oo = atoms.get_distance(0, 3)
vector_oo = atoms.get_distance(0, 3, mic=True, vector=True)
direction_oo = vector_oo / distance_oo
distances = np.linspace(2.0, 6.0, 10)
displacements = distances - distance_oo
potentialenergies = []

for d, dist in zip(displacements, distances):
    atoms_distanced = atoms.copy()
    atoms_distanced[3].position += direction_oo * d
    atoms_distanced[4].position += direction_oo * d
    atoms_distanced[5].position += direction_oo * d
    atoms_distanced.center(vacuum=50.0)
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
