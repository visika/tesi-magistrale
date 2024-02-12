#!/usr/bin/env python3
from ase.io import write
from ase import Atoms
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from mace.calculators import mace_mp

atoms = Atoms("HOH", positions=[(0, 0, 0), (0, 1, 0), (0, 0, 1)])
atoms.center(vacuum=3.0)
calc = mace_mp(model="medium", dispersion=False, default_dtype="float64", device="cpu")
atoms.calc = calc

opt = BFGS(atoms, logfile="optimization.log", trajectory="optimization.traj")
opt.run(fmax=0.0001, steps=1000)
write("final.pdb", atoms)

potentialenergy = atoms.get_potential_energy()

f = open("e_gas.txt", "w")
f.write(str(potentialenergy))
f.close()

for d in [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]:
    vib = Vibrations(atoms, delta=d, name=f"vib_delta={d}")
    vib.run()
    vib.summary(log=f"H2O_delta={d}_summary.txt")
