#!/usr/bin/env python3
import os
import pandas as pd
from ase.io import read
from mace.calculators import mace_mp

structures = [
    "Ih",
    "II",
    "III",
    "IV",
    "VI",
    "VII",
    "VIII",
    "IX",
    "XI",
    "XIII",
    "XIV",
    "XV",
    "XVII",
]

crystal_energies = []

calc = mace_mp(model="medium", dispersion=True, default_dtype="float64", device="cpu")

for s in structures:
    os.chdir(s)
    atoms = read("POSCAR")
    atoms.calc = calc
    nmol = atoms.get_global_number_of_atoms() / 3.0
    energy = atoms.get_potential_energy()
    e_crys = energy / nmol
    crystal_energies.append(e_crys)
    os.chdir("..")

dict = {"structure": structures, "e_crys": crystal_energies}
df = pd.DataFrame(data=dict)
df.to_csv("crystal_energies.csv", index=False)
