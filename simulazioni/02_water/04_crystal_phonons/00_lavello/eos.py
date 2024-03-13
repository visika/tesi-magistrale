# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 23-11-2023
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html

from __future__ import annotations

import numpy as np

from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.eos import EquationOfState
from ase.io import read
from ase.units import kJ
from ase import io
from ase.filters import ExpCellFilter, FrechetCellFilter
from ase.optimize import LBFGS
from mlip_calculators import choose_calculator


import argparse
import pathlib

cli = argparse.ArgumentParser()
cli.add_argument(
    "--device",
    type=str,
    default="cuda",
    choices=["cuda", "cpu", "mps"],
    help=" provide device to run. default %(default)s",
)

cli.add_argument(
    "--precision",
    type=str,
    default="float64",
    choices=["float64", "float32"],
    help=" provide precision to run. default %(default)s",
)

cli.add_argument("--d3", action="store_true", help="enable d3? default %(default)s")

cli.add_argument(
    "--system",
    type=str,
    default=None,
    help="system coordinates to simulate, ase readable, default %(default)s",
)

cli.add_argument(
    "--min_scale",
    type=float,
    default=0.95,
    help="minimum scale in 0-1, default %(default)s",
)

cli.add_argument(
    "--max_scale", type=float, default=1.05, help="max scale in 1-, default %(default)s"
)

cli.add_argument(
    "--images",
    type=int,
    default=10,
    help="minimum scale in integers, default %(default)s",
)
cli.add_argument(
    "--model",
    type=str,
    default="",
    help="provide the mace model file. default: %(default)s",
)

cli.add_argument(
    "--minimize_fmax",
    type=float,
    default=0.01,
    help="provide minimizer precision for forces in eV/A, default %(default)s",
)


def optimize(a: Atoms, fmax: float, optf: str) -> None:
    ua = FrechetCellFilter(a, hydrostatic_strain=True)
    dyn = LBFGS(ua)
    dyn.run(fmax=fmax)
    io.write(optf, a, write_info=True, format="extxyz")
    io.write("opt.cif", a, format="cif")


args = cli.parse_args()

model = args.model
device = args.device
system = args.system
precision = args.precision
d3 = args.d3
n = args.images
smin = args.min_scale
smax = args.max_scale
fmax = args.minimize_fmax

sys = read(system)
sysname = pathlib.Path(system).stem
calculator = choose_calculator(
    architecture="mace_mp",
    model=model,
    dispersion=d3,
    default_dtype=precision,
    device=device,
)
sys.calc = calculator

m = model
if d3:
    m = f"{model}-d3"
optimize(sys, fmax, f"{sysname}-{m}-opt.xyz")

cell = sys.get_cell()

traj = Trajectory(f"{sysname}.traj", "w")
for x in np.linspace(smin, smax, n):
    sys.set_cell(cell * x, scale_atoms=True)
    sys.get_potential_energy()
    traj.write(sys)

configs = read(f"{sysname}.traj@0:{n}")  # read 5 configurations
# Extract volumes and energies:
volumes = [sys.get_volume() for sys in configs]
energies = [sys.get_potential_energy() for sys in configs]
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
print(B / kJ * 1.0e24, "GPa")
eos.plot(f"{sysname}-{m}-eos.png")
