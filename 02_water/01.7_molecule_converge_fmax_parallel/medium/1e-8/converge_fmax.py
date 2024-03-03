#!/usr/bin/env python3
import argparse
from ase import Atoms
from ase.io import write
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from mace.calculators import mace_mp

parser = argparse.ArgumentParser()
parser.add_argument("--dispersion", action="store_true")
parser.add_argument("--model")
parser.add_argument("--fmax")
args = parser.parse_args()

print("fmax: {float(args.fmax)}")
print(f"model: {args.model}")
print(f"dispersion: {args.dispersion}")

# atoms = ase.build.molecule("H2O")
atoms = Atoms("HOH", positions=[(0, 0, 0), (0, 1, 0), (0, 0, 1)])
atoms.calc = mace_mp(
    model=args.model,
    dispersion=args.dispersion,
    default_dtype="float64",
    device="cpu",
)
opt = BFGS(atoms, logfile="optimization.log", trajectory="optimization.traj")
opt.run(fmax=float(args.fmax), steps=1000)
write("final.xyz", atoms)
potentialenergy = atoms.get_potential_energy()
f = open("e_gas.txt", "w")
f.write(str(potentialenergy))
f.close()
for d in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]:
    vib = Vibrations(atoms, delta=d, nfree=4, name=f"vib_delta={d}")
    vib.run()
    vib.summary(log=f"H2O_delta={d}_summary.txt")
    vib.clean()
