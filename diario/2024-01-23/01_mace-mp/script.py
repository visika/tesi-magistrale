#!/usr/bin/env python3
from mace.calculators import mace_mp
from ase import build

atoms = build.molecule("H2O")
calc = mace_mp(model="medium", dispersion=False, default_dtype="float32", device="cuda")
atoms.calc = calc
print(atoms.get_potential_energy())
