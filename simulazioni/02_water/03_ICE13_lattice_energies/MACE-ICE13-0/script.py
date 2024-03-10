#!/usr/bin/env python3
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

calc = mace_mp(
    model="/ibiscostorage/VirtualMatterLab/MACE-ICE13/MACE-ICE13.model",
    dispersion=True,
    default_dtype="float64",
    device="cuda",
)

for s in structures:
    atoms = read(f"/ibiscostorage/mmollo/tesi-magistrale/strutture/ICE13/{s}/POSCAR")
    atoms.calc = calc
    nmol = atoms.get_global_number_of_atoms() / 3.0
    energy = atoms.get_potential_energy()
    e_crys = energy / nmol
    crystal_energies.append(e_crys)

dict = {"structure": structures, "e_crys": crystal_energies}
df = pd.DataFrame(data=dict)
df.to_csv("crystal_energies.csv", index=False)
