#!/usr/bin/env python3
from ase.io import read
from ase.vibrations import Vibrations
from mace.calculators.mace import MACECalculator


calculator = MACECalculator(
    "/ibiscostorage/VirtualMatterLab/MACE-ICE13/MACE-ICE13-1.model",
    default_dtype="float64",
    device="cpu",
)
atoms = read(
    "/ibiscostorage/mmollo/tesi-magistrale/strutture/patridge1997_gas_phase/POSCAR"
)
atoms.calc = calculator
potentialenergy = atoms.get_potential_energy()
f = open("e_gas.txt", "w")
f.write(str(potentialenergy))
f.close()

for delta in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]:
    vib = Vibrations(atoms, delta=delta, nfree=4, name=f"vib_delta={delta}")
    vib.run()
    vib.summary()
    vib.summary(log=f"H2O_delta={delta}_summary.txt")
    vib.clean()
