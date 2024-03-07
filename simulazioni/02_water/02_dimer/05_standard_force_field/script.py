#!/usr/bin/env python3
from ase.io import read, write
from ase.optimize import BFGS
from ase.calculators.harmonic import HarmonicForceField, HarmonicCalculator

ref_atoms = read("init.xyz")

hff = HarmonicForceField(ref_atoms=ref_atoms)
