#!/usr/bin/env python3
from ase.vibrations.vibrations import Vibrations
from mace.calculators import mace_mp
from ase.io import read

atoms = read("relax-monomer-hard/final.pdb")
# atoms.center(vacuum=3.0)
calc = mace_mp(model="large", dispersion=False, default_dtype="float64", device="cuda")
atoms.calc = calc

vib = Vibrations(atoms)
vib.run()
vib.summary(log="H2O_MACE_summary.txt")
vib.write_mode()
