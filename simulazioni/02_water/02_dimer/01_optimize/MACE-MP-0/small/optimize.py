#!/usr/bin/env python3
import argparse

from ase.io import read, write
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from mace.calculators import mace_mp

parser = argparse.ArgumentParser()
parser.add_argument("--model")
args = parser.parse_args()

print(f"model: {args.model}")

atoms = read("../init.xyz")
atoms.calc = mace_mp(
    model=args.model, dispersion=True, default_dtype="float64", device="cuda"
)

opt = BFGS(atoms, logfile="optimization.log", trajectory="optimization.traj")
opt.run(fmax=1e-8, steps=1000)
write("final.xyz", atoms)

for delta in [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]:
    vib = Vibrations(atoms, delta=delta, nfree=4, name=f"vib_delta={delta}")
    vib.run()
    vib.summary(log=f"dimer_delta={delta}_summary.txt")
    vib.clean()
