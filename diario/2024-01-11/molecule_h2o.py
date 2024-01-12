#!/usr/bin/env python3
from ase.build import molecule
from ase.visualize import view

atoms = molecule("H2O", vacuum=3.0)
view(atoms)
