# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 23-11-2023
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html

from ase import io
from ase.mep import DyNEB, NEB
from ase.optimize import LBFGS
from ase.mep import NEBTools

from mace.calculators import mace_mp
import argparse
from ase.constraints import FixAtoms
import pathlib

cli=argparse.ArgumentParser()
cli.add_argument(
  "--d3",
  action="store_true",
  help = 'enable d3? default %(default)s'
)

cli.add_argument(
  "--model",
  type=str,
  default='',
  help = 'provide the mace model file. default: %(default)s'
)

cli.add_argument(
  "--precision",
  type=str,
  default='float64',
  choices = ['float64','float32'],
  help = ' provide precision to run. default %(default)s'
)

cli.add_argument(
  "--method",
  type=str,
  default='dyneb',
  choices = ['dyneb','neb'],
  help = 'method . default %(default)s'
)

cli.add_argument(
  "--A",
  type=str,
  default='stateA.xyz',
  help = 'coordinates for state A. default: %(default)s'
)

cli.add_argument(
  "--B",
  type=str,
  default='stateB.xyz',
  help = 'coordinates for state B. default: %(default)s'
)

cli.add_argument(
  "--traj",
  type=str,
  default='traj',
  help = 'trajectory for ne. default: %(default)s'
)

cli.add_argument(
  "--device",
  type=str,
  default='cuda',
  choices = ['cuda','cpu','mps'],
  help = ' provide device to run. default %(default)s'
)

cli.add_argument(
  "--nimages",
  type=int,
  default=15,
  help = 'number of images. default: %(default)s'
)

cli.add_argument(
  "--fmax",
  type=float,
  default=0.01,
  help = 'fmax for optimisation. default: %(default)s'
)

cli.add_argument(
  "--fixed_atoms",
  nargs="*",
  type=int,
  default=[],
  help = ' list of atoms to fix... default %(default)s'
)


args = cli.parse_args()
precision = args.precision
device = args.device
with_d3 = args.d3
fmax = args.fmax
nimages = args.nimages
stateA = args.A
stateB = args.B
fixed_atoms = args.fixed_atoms
model = args.model
traj = args.traj
method = args.method


m = model.split("/")[-1].replace(".model","")
if with_d3:
   m += "-d3"
m += f"-{method}"
sA = pathlib.Path(stateA).stem
sB = pathlib.Path(stateB).stem
print(f"enable d3: {with_d3}")
print(f"model: {model}")
print(f"method: {method}")
initial = io.read(stateA)
if len(fixed_atoms) > 0:
  c = FixAtoms(indices=fixed_atoms)
  initial.set_constraint(c)

initial.calc = mace_mp(model=model, device=device, default_dtype=precision, dispersion=with_d3)

g_in = LBFGS(initial, trajectory=f"{sA}-{m}.traj")
g_in.run(fmax=fmax)
io.write(f"{sA}-opt-{m}.cif",initial)
final = io.read(stateB)

if len(fixed_atoms) > 0:
  c = FixAtoms(indices=fixed_atoms)
  final.set_constraint(c)
final.calc = mace_mp(model=model, device=device, default_dtype=precision, dispersion=with_d3)
g_fin = LBFGS(final, trajectory=f"{sB}-{m}.traj")
g_fin.run(fmax=fmax)
io.write(f"{sB}-opt-{m}.cif",final)

initial = io.read(f"{sA}-{m}.traj")
final = io.read(f"{sB}-{m}.traj")
io.write(f"{sA}-opt-{m}.xyz",initial,format="extxyz")
io.write(f"{sB}-opt-{m}.xyz",final,format="extxyz")

images = [initial]
images += [initial.copy() for i in range(nimages)]
images += [final]
match method:
  case 'neb':
    neb = NEB(images, k=0.1, climb=True,method='string')
  case 'dyneb':
    neb = DyNEB(images, fmax=fmax, dynamic_relaxation=True, climb=True, scale_fmax=1.2)

neb.interpolate(method='idpp')

for i in range(len(images)):
   images[i].calc = mace_mp(model=model, device=device, default_dtype=precision, dispersion=with_d3)

optimizer = LBFGS(neb, trajectory=f"{traj}-{m}.traj")
optimizer.run(fmax=fmax)
io.write(f"{traj}-{m}.xyz",images, write_info=False)
images = io.read(f"{traj}-{m}.traj@-{nimages+2}:")
nebtools = NEBTools(images[1:-1])
Ef, dE = nebtools.get_barrier()

print(Ef,dE)
max_force = nebtools.get_fmax()

fig = nebtools.plot_band()
fig.savefig(f'{traj}-{m}.png')
