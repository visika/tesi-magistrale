#!/usr/bin/env python3
from ase import Atoms
from ase.io import write
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from mace.calculators.mace import MACECalculator


calculator = MACECalculator(
    "/ibiscostorage/VirtualMatterLab/MACE-ICE13/MACE-ICE13-1.model",
    default_dtype="float64",
    device="cuda",
)
atoms = Atoms(
    symbols="H2O",
    positions=[(0, 0, 0), (0, 1, 0), (0, 0, 1)],
    pbc=False,
    calculator=calculator,
)
opt = BFGS(atoms, logfile="optimization.log", trajectory="optimization.traj")
if opt.run(fmax=1e-8, steps=1000):
    print("Converged!")
else:
    print("Not converged!")
write("final.xyz", images=atoms)
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
