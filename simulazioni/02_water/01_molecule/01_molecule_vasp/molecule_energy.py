#!/usr/bin/env python3
from ase.io import read
from ase.calculators.vasp import Vasp

calc = Vasp(txt="vasp.out", ismear=0, sigma=0.01)

atoms = read("final.xyz")
atoms.center(vacuum=10)
atoms.pbc = True
atoms.calc = calc

potentialenergy = atoms.get_potential_energy()

f = open("e_gas.txt", "w")
f.write(str(potentialenergy))
f.close()
