#!/usr/bin/env python3
from ase.io import read
from ase.phonons import Phonons
from mace.calculators import mace_mp

atoms = read("relax-full/final_1.0GPa.pdb")
calc = mace_mp(model="medium", dispersion=False, default_dtype="float64", device="cuda")
atoms.calc = calc

# Supercell number
N = 2
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.01)
ph.run()
# Find the json data in the phonon folder
