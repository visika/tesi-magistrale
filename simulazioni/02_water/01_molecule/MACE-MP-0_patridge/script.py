#!/usr/bin/env python3
from ase.io import read
from mace.calculators import mace_mp

atoms = read(
    "/ibiscostorage/mmollo/tesi-magistrale/strutture/patridge1997_gas_phase/POSCAR"
)

calculator = mace_mp(
    model="medium", dispersion=True, default_dtype="float64", device="cuda"
)

atoms.calc = calculator

potentialenergy = atoms.get_potential_energy()

with open("e_gas.txt", "w") as f:
    f.write(str(potentialenergy))
