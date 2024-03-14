#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 23-11-2023
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html

from ase import units
from ase.md.npt import NPT
from ase.md.langevin import Langevin
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution,Stationary,ZeroRotation
from ase.io import read, write
from ase.optimize import LBFGS,MDMin
from ase.filters import FrechetCellFilter

import numpy as np
from mace.calculators import mace_mp
import datetime as clock
import sys
import os

from ase.constraints import FixAtoms, FixBondLengths
import argparse
import pathlib
import shutil

densfact = (units.m/1.0e2)**3/units.mol

def optimize(a,fmax,optimize_cell,output):
  x = os.getpid()
  traj = f"tmp{x}.traj"
  if optimize_cell:
    ua = FrechetCellFilter(a,hydrostatic_strain=True)
    dyn = LBFGS(ua,trajectory=traj,logfile=opt_log)
  else:
    dyn = LBFGS(a,trajectory=traj,logfile=opt_log)

  print(f"optimization started at {clock.datetime.now()}")
  dyn.run(fmax=fmax)
  write(output,a,format="extxyz")
  x = read(traj,index=":")
  write(opt_traj,x,format="extxyz")
  write(pathlib.PurePosixPath(output).stem+'.cif',a,format="cif")
  try:
     p = pathlib.Path(traj)
     p.unlink()
  except:
     pass

  print(f"optimization  finished at {clock.datetime.now()}")

def rotate_restart(files,restarts_to_keep):

  if len(files) > restarts_to_keep:
     try:
       p = pathlib.Path(files.pop(0))
       p.unlink()
     except:
       pass

cli=argparse.ArgumentParser()

cli.add_argument(
  "--barostat_time",
  type=float,
  default=75.0,
  help = 'provide barostat, in fs, default %(default)s'
)

cli.add_argument(
  "--bulk_modulus",
  type=float,
  default=2.0,
  help = 'provide bulk modulus for npt, in GPa, default %(default)s'
)

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
   "--ensemble",
   type=str,
   default = 'nve',
   choices = ['nve','npt','nvt','nvt-nh'],
   help = "choose the ensemble you want to simulate. default %(default)s"
)

cli.add_argument(
  "--equil_steps",
  type=int,
  default=0,
  help = 'number of equilibration steps in simulation default %(default)s'
)


cli.add_argument(
  "--fixed_atoms",
  nargs="*",
  type=int,
  default=[],
  help = ' list of atoms to fix... default %(default)s'
)

cli.add_argument(
  "--fixed_bond",
  nargs="+",
  type=int,
  default=[],
  action = "append",
  help = ' list of bonds to fix, add as many --fixed_bond as needed default %(default)s'
)

cli.add_argument(
  "--fully_flexible",
  action="store_true",
  help = 'controls mask for NPT simulations? default %(default)s only lengths of cell change'
)

cli.add_argument(
  "--group_atoms",
  nargs="*",
  type=int,
  default=[],
  help = ' list of atoms to print more often, hardcoded to 10 steps now... default %(default)s'
)

cli.add_argument(
  "--group_every",
  type=int,
  default=100,
  help = 'frequency to print a group of atoms, default %(default)s'
)

cli.add_argument(
  "--group_traj",
  type=str,
  default='group',
  help = ' group of atoms traj name base. default %(default)s'
)


cli.add_argument(
  "--lang_friction",
  type=float,
  default=0.005,
  help = 'provide thermfstat, in fs^-1, default %(default)s'
)

cli.add_argument(
  "--log",
  type=str,
  default='',
  help = ' provide log file name. default %(default)s'
)

cli.add_argument(
  "--minimize",
  action="store_true",
  help = 'shall I minimize during equilibration? default %(default)s'
)

cli.add_argument(
  "--minimize_every",
  type=int,
  default=-1,
  help = 'frequency to minimize every nth step during equilibration, default %(default)s'
)

cli.add_argument(
  "--minimize_fmax",
  type=float,
  default=0.01,
  help = 'maximum force on atoms for optimisation, in eV/A, default %(default)s'
)

cli.add_argument(
  "--model",
  type=str,
  default='',
  help = 'provide the mace model file. default: %(default)s'
)

cli.add_argument(
  "--output_every",
  type=int,
  default=100,
  help = 'frequency to output log, default %(default)s'
)

cli.add_argument(
  "--optimize_cell",
  action="store_true",
  help = 'shall I optimize cell in addition to positions? default %(default)s'
)

cli.add_argument(
  "--opt_struct",
  type=str,
  default='',
  help = 'write optimised structure to. default: %(default)s'
)

cli.add_argument(
  "--opt_traj",
  type=str,
  default='',
  help = 'write opt trajectory. default: %(default)s'
)

cli.add_argument(
  "--opt_log",
  type=str,
  default='',
  help = 'write optimization log. default: %(default)s'
)

cli.add_argument(
  "--pressure",
  type=float,
  default=0.0,
  help = 'pressure in bar, default %(default)s'
)

cli.add_argument(
  "--rescale_every",
  type=int,
  default=13,
  help = 'frequency to rescale velocities to target, default %(default)s'
)

cli.add_argument(
  "--rescale_velocities",
  action="store_true",
  help = 'rescale velocities? default %(default)s'
)

cli.add_argument(
  "--restart",
  action="store_true",
  help = 'is the simulation a restart? default %(default)s'
)

cli.add_argument(
  "--restart_auto",
  action="store_true",
  help = 'restarts a simulation auto from were was left  default %(default)s'
)

cli.add_argument(
  "--restart_rotate",
  action="store_true",
  help = 'shall we rotate the restart files? default %(default)s'
)

cli.add_argument(
  "--restart_every",
  type=int,
  default=1000,
  help = 'restart info every x steps:  default %(default)s'
)

cli.add_argument(
  "--restart_keep",
  type=int,
  default=4,
  help = 'restart structures to keep:  default %(default)s'
)

cli.add_argument(
  "--seed",
  type=int,
  default=2023,
  help = 'random seed for langevin, default %(default)s'
)

cli.add_argument(
  "--steps",
  type=int,
  default=0,
  help = 'number of steps in simulation, default %(default)s'
)

cli.add_argument(
  "--system",
  type=str,
  default=None,
  help = 'system coordinates to simulate, ase readable, default %(default)s'
)

cli.add_argument(
  "--system_name",
  type=str,
  default='',
  help = 'name of the system to simulate, determines defaults of outputs if not overwritten, default %(default)s'
)

cli.add_argument(
  "--traj",
  type=str,
  default='',
)

cli.add_argument(
  "--traj_append",
  action="store_true",
  help = 'append to the trajectory? default %(default)s'
)

cli.add_argument(
  "--traj_start",
  type=int,
  default=0,
  help = 'start writing trajectory at step:  default %(default)s'
)

cli.add_argument(
  "--traj_every",
  type=int,
  default=100,
  help = 'trajectory write every:  default %(default)s'
)

cli.add_argument(
  "--thermostat_time",
  type=float,
  default=50.0,
  help = 'provide thermostat, in fs, default %(default)s'
)

cli.add_argument(
  "--timestep",
  type=float,
  default=1.0,
  help = 'provide timestep for integrator, in fs, default %(default)s'
)

cli.add_argument(
  "--temp_start",
  type=float,
  default=-1,
  help = 'provide temperature to start heating, negative means disabled, in K, default %(default)s'
)

cli.add_argument(
  "--temp_end",
  type=float,
  default=-1,
  help = 'provide temperature to end heating, negative means disabled, in K, default %(default)s'
)

cli.add_argument(
  "--temp_step",
  type=float,
  default=-1,
  help = 'step to increase the Temperature duting heating, in K, default %(default)s'
)

cli.add_argument(
  "--temp_time",
  type=float,
  default=-1,
  help = 'time between heating steps, in fs, default %(default)s'
)

cli.add_argument(
  "--temperature",
  type=float,
  default=300,
  help = 'provide Temperature in K, default %(default)s'
)


args = cli.parse_args()

model = args.model
device = args.device
system = args.system
precision = args.precision
sysname = args.system_name
if sysname == '':
   sysname = pathlib.Path(system).stem

timestep = args.timestep
steps = args.steps
equil_steps = args.equil_steps
T = args.temperature

rescale_every = args.rescale_every
rescale_velocities = args.rescale_velocities

output_every = args.output_every

fixed_atoms = args.fixed_atoms
fixed_bonds = args.fixed_bond

traj = args.traj
traj_append = args.traj_append
traj_start = args.traj_start
traj_every = args.traj_every

with_d3 = args.d3

is_restart = args.restart
restart_every = args.restart_every
restarts_to_keep = args.restart_keep
restart_rotate = args.restart_rotate
restart_auto = args.restart_auto

ensemble = args.ensemble
ttime = args.thermostat_time
ptime = args.barostat_time
bulk_modulus = args.bulk_modulus
ext_pressure = args.pressure
lang_friction = args.lang_friction

log_file = args.log

minimize = args.minimize
minimize_every = args.minimize_every
fmax = args.minimize_fmax
optimize_cell = args.optimize_cell
optimized_struct = args.opt_struct
opt_log = args.opt_log
opt_traj = args.opt_traj

group=args.group_atoms
group_every = args.group_every
group_traj= args.group_traj

temp_start = args.temp_start
temp_end = args.temp_end
temp_step = args.temp_step
temp_time = args.temp_time

fully_flexible = args.fully_flexible
mask = np.eye(3)
if fully_flexible:
   mask = np.ones((3,3))
seed = args.seed
if restart_auto:
   traj_append = True
   is_restart = True


if log_file == '':
  log_file =f"{sysname}-{ensemble}-{T}-md.log"
if traj == '':
  traj =f"{sysname}-{ensemble}-{T}-traj.xyz"
out_file = f"{sysname}-{ensemble}-{T}-md.out"

of = open(out_file,"w")
print("#"+" ".join(sys.argv),file=of)
print(f"{sysname} model setup started at {clock.datetime.now()}",file=of)

if model =='':
   print("you need to provide a model! use --model name or path")
   print("you need to provide a model! use --model name or path",file=of)
   exit()
m = pathlib.Path(model).stem
if with_d3:
    m = f"{model}-d3"
if opt_log == '':
  opt_log =f"{sysname}-opt.log"
if optimized_struct == '':
  optimized_struct = f"{sysname}-{m}-opt.xyz"
if opt_traj == '':
  opt_traj = f"{sysname}-opt-traj.xyz"
if len(group)>0:
  group_traj = f"{sysname}-{group_traj}-traj.xyz"

restart_stem = f"res-{sysname}-{ensemble}"

offset = 0
# try to guess the restart
if is_restart:
  try:
    with open(log_file,'r') as f:
      # we addume logs are small
      last_line = f.readlines()[-1].split()
      try:
        offset = int(last_line[0])
      except:
        pass
  except:
    pass

print_header = False
if restart_auto:
   try:
     possible_restart = pathlib.Path(".").glob(f"{restart_stem}*.xyz")
     last_restart = list(sorted(possible_restart,key=os.path.getmtime))[-1]
     shutil.copy(last_restart,pathlib.Path(f"{sysname}-{last_restart}"))
     system = last_restart
   except:
      print_header = True
      pass

log=open(log_file,"a" if is_restart else "w")


mace_calc = mace_mp(model=model, device=device, default_dtype=precision, dispersion=with_d3)
asystem = read(system)
if "masses" not in asystem.arrays.keys():
    asystem.set_masses()
n = len(asystem)
for k in ['spacegroup', 'unit_cell', 'occupancy']:
  if k in asystem.info.keys():
     del asystem.info[k]

asystem.info["real_time"] = clock.datetime.now()
if temp_start > 0:
   T = temp_start

print(f"""
system: {system}
No of atoms: {n}
model: {model}
d3 enabled: {with_d3}
equilibration steps: {equil_steps}
output every: {output_every}
timestep: {timestep} fs
log file name: {log_file}
fixed atoms: {fixed_atoms}
fixed bonds: {fixed_bonds}
device: {device}
restart: {is_restart}
seed: {seed}
""", file=of)
if is_restart:
  print(f"restart_every: {restart_every}", file=of)

print(f"""
trajectory name: {traj}
trajectory start:{traj_start}
trajectory every: {traj_every}
trajectory append: {traj_append}
""", file=of)

if temp_start > 0:
   steps = int(temp_time/timestep)
   print(f"""Temperature: ramping
temperature start: {temp_start} K
temperature end: {temp_end} K
temperature increase: {temp_step} K
temperature time perios: {temp_time} fs
steps: {steps} steps per T
""", file=of)
else:
   print(f"""Temperature: {T}
steps: {steps}
""", file=of)

if minimize:
  print(f"""
minimize: {minimize}
minimize every: {minimize_every} steps
minimise fmax: {fmax}
optimize cell: {optimize_cell}
optimized structure: {optimized_struct}
""",file=of)

print(f"ensemble: {ensemble}", file=of)
match ensemble:
  case  "npt":
   print(f"""
thermostat time: {ttime} fs
barostat time: {ptime} fs
bulk modules: {bulk_modulus} GPa
pressure: {ext_pressure} bar
""", file=of)
  case "nvt"|'nvt-nh':
    print(f"friction: {lang_friction}",file=of)

print(f"rescale velocities: {rescale_velocities}",file=of)
if rescale_velocities:
    print(f"rescale every: {rescale_every}  steps",file=of)
print(f"offset: {offset}",file=of)
log.flush()
of.flush()

if len(fixed_atoms) > 0:
   c = FixAtoms(indices=fixed_atoms)
   asystem.set_constraint(c)
   print("Atoms fixed:",file=of)
   for i,a in enumerate([ asystem[i] for i in c.get_indices()]):
     print(f"  - {i}: {a.symbol} {a.index}",file=of)

if len(fixed_bonds) > 0:
   c = FixBondLengths(fixed_bonds)
   asystem.set_constraint(c)
   bonds = c.get_indices()
   print("Bonds fixed:",file=of)
   for i,b in   enumerate([ (asystem[bonds[i]], asystem[bonds[i+1]]) for i in range(0,len(bonds),2)]):
     print(f"  - {i}: {b[0].symbol}({b[0].index}) - {b[1].symbol}({b[1].index})",file=of)


asystem.calc = mace_calc
if len(group) == 1:
    k = group[0]
    if group[0]< 0:
        group = [ i for i in range(len(asystem)+k, len(asystem))]
    print(f"dumping every {group_every} in {group_traj} the group of atoms: {group}")

print(f"{sysname}: model setup finished at {clock.datetime.now()}",file=of)

np.random.seed(seed)
match ensemble:
  case "npt":
    ttime *= units.fs
    ptime *= units.fs
    BulkModulus=bulk_modulus*units.GPa
    pfactor=ptime**2*BulkModulus
    external_stress = ext_pressure * units.bar
  case "nvt":
    lang_friction /= units.fs
  case "nvt-nh":
    ttime *= units.fs

def md(T: float, p: float, restart: bool, traj_append: bool, offset: int) -> None:
  match ensemble:
    case "npt":
      dyn = NPT(asystem, timestep=timestep*units.fs, temperature_K=T,ttime=ttime,pfactor=pfactor,
	mask=mask, append_trajectory=traj_append,
	externalstress=external_stress)
    case "nvt-nh":
      dyn = NPT(asystem, timestep=timestep*units.fs, temperature_K=T,ttime=ttime,pfactor=None, append_trajectory=traj_append)
    case "nvt":
      dyn = Langevin(asystem, timestep=timestep*units.fs, temperature_K=T, friction=lang_friction, append_trajectory=traj_append)
    case "nve":
      dyn = VelocityVerlet(asystem, timestep=timestep*units.fs, append_trajectory=traj_append)
  restart_files =[]
  def write_frame():
    if dyn.nsteps > traj_start and dyn.nsteps % traj_every == 0:
      dyn.atoms.write(traj, write_info=True, columns=['symbols','positions','momenta','masses'], append=True)
    epot = dyn.atoms.get_potential_energy() / n
    ekin = dyn.atoms.get_kinetic_energy() / n
    cT = ekin / (1.5 * units.kB)
    time = offset*timestep + dyn.get_time()/units.fs
    step = offset + dyn.nsteps
    dyn.atoms.info['time_fs'] = time
    dyn.atoms.info['step'] = step
    t = clock.datetime.now()
    real_time = t - dyn.atoms.info['real_time']
    dyn.atoms.info['real_time'] = t
    try:
      density = np.sum(dyn.atoms.get_masses())/dyn.atoms.get_volume()*densfact
      dyn.atoms.info['density'] = density
      volume = dyn.atoms.get_volume()
      pressure = -np.trace(dyn.atoms.get_stress(include_ideal_gas=True, voigt=False))/3/units.bar
      pressure_tensor = - dyn.atoms.get_stress(include_ideal_gas=True, voigt=True)/units.bar
    except:
      volume = 0
      pressure = 0
      density = 0
      pressure_tensor = np.zeros(6)
    print_stat=f"{step:10d} {real_time.total_seconds():.3f} {time:13.2f} {epot:.3e} {ekin:.3e} {cT:.3f} {epot+ekin:.3e} {density:.3f} {volume:.3e} {pressure:.3e} {pressure_tensor[0]:.3e} {pressure_tensor[1]:.3e} {pressure_tensor[2]:.3e} {pressure_tensor[3]:.3e} {pressure_tensor[4]:.3e} {pressure_tensor[5]:.3e}"

    match ensemble:
      case "npt":
        print_stat += f" {ext_pressure} {T}"
      case "nvt"| "nvt-nh":
        print_stat += f" {T}"
    print(print_stat,file=log)
    if step%restart_every == 0:
      ff = f"{restart_stem}-{T:.2f}K-{step}.xyz"
      write(ff,asystem,write_info=True, columns=['symbols','positions','momenta','masses'])
      if restart_rotate:
        restart_files.append(ff)
        rotate_restart(restart_files,restarts_to_keep)
      log.flush()

  def write_group():
      write(group_traj, dyn.atoms[group],write_info=True, columns=['symbols','positions','momenta','masses'], append=True)

  def reset_velocities():
    if dyn.nsteps < equil_steps:
      MaxwellBoltzmannDistribution(asystem, temperature_K=T)
      Stationary(asystem)
      ZeroRotation(asystem)

  def optimize_structure():
    if dyn.nsteps < equil_steps:
       optimize(asystem,fmax,optimize_cell,optimized_struct)


  if not restart or print_header:
    print_h = '    Step  | real time[s] |     Time [fs] |   Epot/N [eV] | Ekin/N [eV] |  T [K] | Etot/N [eV] | Density [g/cm^3] |Volume [A^3] | Pressure [bar] |   Pxx [bar] |   Pyy [bar] |   Pzz[bar] |   Pyz[bar] |   Pxz[bar] |   Pxy[bar]'
    match ensemble:
      case "nvt"|"nvt-nh":
        print_h += ' |Target T [K]'
      case "npt":
        print_h += ' |Target Pressure[bar] |  T [K]'
    print(print_h, file=log)

  dyn.attach(write_frame, interval=output_every)
  if rescale_velocities:
    dyn.attach(reset_velocities, interval=rescale_every)

  if minimize and  minimize_every > 0:
    dyn.attach(optimize_structure, interval=minimize_every)
  if len(group)>0:
    dyn.attach(write_group, interval=group_every)

  dyn.run(steps)

if minimize and not is_restart:
  optimize(asystem,fmax,optimize_cell,optimized_struct)
  MaxwellBoltzmannDistribution(asystem, temperature_K=T)
  Stationary(asystem)
  ZeroRotation(asystem)

if not minimize and not is_restart:
  MaxwellBoltzmannDistribution(asystem, temperature_K=T)
  Stationary(asystem)
  ZeroRotation(asystem)

restart = is_restart

t1 = clock.datetime.now()
print(f"{sysname} start running MD at: {t1}",file=of)
if temp_start > 0:
  while T <= temp_end:
    a = clock.datetime.now()
    print(f"{sysname} MD for {T=} started at {a}",file=of)
    p = external_stress
    md(T,p,restart,traj_append,offset)
    b = clock.datetime.now()
    print(f"{sysname} MD for {T=} ended at {b}, elpased time {b-a}",file=of)
    of.flush()
    T += temp_step
    restart = True
    traj_append = True
    offset += steps
else:
    if steps > 0 :
      p = ext_pressure
      md(T,p,restart,traj_append,offset)
    else:
      print(f"Potential energy of the systems is: {asystem.get_potential_energy()}")

t2 = clock.datetime.now()
print(f"{sysname}: MD finished at {t2}, elapsed time {t2-t1}",file=of)
of.close()
