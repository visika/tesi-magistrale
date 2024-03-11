# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 31-01-2024
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html

from __future__ import annotations

import numpy as np

from ase import Atoms
from ase.io import read
from ase import io
from ase.filters import FrechetCellFilter
from ase.optimize.precon import PreconFIRE,Exp
from ase.optimize import LBFGS
from mlip_calculators import choose_calculator


import argparse
import pathlib

cli=argparse.ArgumentParser()
cli.add_argument(
  "--device",
  type=str,
  default='cuda',
  choices = ['cuda','cpu','mps'],
  help = ' provide device to run. default %(default)s'
)

cli.add_argument(
  "--precision",
  type=str,
  default='float64',
  choices = ['float64','float32'],
  help = ' provide precision to run. default %(default)s'
)

cli.add_argument(
  "--d3",
  action="store_true",
  help = 'enable d3? default %(default)s'
)

cli.add_argument(
  "--system",
  type=str,
  default=None,
  help = 'system coordinates to simulate, ase readable, default %(default)s'
)

cli.add_argument(
  "--optimize_cell",
  action="store_true",
  help = 'enable cell optimization? default %(default)s'
)

cli.add_argument(
  "--model",
  type=str,
  default='',
  help = 'provide the mace model file. default: %(default)s'
)

cli.add_argument(
  "--minimize_fmax",
  type=float,
  default=0.01,
  help = 'provide minimizer precision for forces in eV/A, default %(default)s'
  )

def optimize(a: Atoms ,fmax : float , optcell: bool, optf: str) -> None:
  #ua = FrechetCellFilter(a,hydrostatic_strain=True)
  #dyn = PreconFIRE(ua,precon=Exp(A=3), use_armijo=True,trajectory="opt.traj")
  dyn = LBFGS(a)
  dyn.run(fmax=fmax)
  io.write(optf,a,write_info=True, format="extxyz")
  io.write('opt.cif',a,format="cif")


args = cli.parse_args()

model = args.model
device = args.device
system = args.system
precision = args.precision
d3 = args.d3
fmax = args.minimize_fmax
optimize_cell = args.optimize_cell

sys = read(system)
sysname = pathlib.Path(system).stem
calculator =choose_calculator(architecture='mace_mp', model=model, dispersion=d3, default_dtype=precision, device=device)
sys.calc = calculator


optimize(sys,fmax,optimize_cell, f"{sysname}-opt.xyz")

print(f"E_config= {sys.get_potential_energy()} eV")
if optimize_cell:
    print(f"cell {sys.get_cell()}")

