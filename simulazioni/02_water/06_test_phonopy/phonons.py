from __future__ import annotations

from pymatgen.core import Structure
import phonopy as phonons
from ase import Atoms
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

from phonopy.structure.atoms import PhonopyAtoms
from ase.io import read, write
from mace.calculators import mace_mp

import matplotlib.pyplot as plt

from ase.optimize import LBFGS, FIRE, MDMin
from ase.optimize.precon import PreconFIRE, Exp

from ase.filters import FrechetCellFilter
import datetime as clock
import pathlib
import os
import argparse

cli = argparse.ArgumentParser()
cli.add_argument("--d3", action="store_true", help="enable d3? default %(default)s")

cli.add_argument(
    "--dos",
    action="store_true",
    help="enable dos/pdos calcualtion? default %(default)s",
)

cli.add_argument(
    "--model",
    type=str,
    default="small",
    help="provide the mace model file. default: %(default)s",
)

cli.add_argument(
    "--precision",
    type=str,
    default="float64",
    choices=["float64", "float32"],
    help=" provide precision to run. default %(default)s",
)

cli.add_argument(
    "--system",
    type=str,
    default="Si.cif",
    help="system coordinates to simulate, ase readable, default %(default)s",
)

cli.add_argument(
    "--minimize_fmax",
    type=float,
    default=0.001,
    help="provide minimizer precision for forces in eV/A, default %(default)s",
)

cli.add_argument(
    "--displacement",
    type=float,
    default=0.01,
    help="provide displacement in A for force constants calculation, default %(default)s",
)

cli.add_argument(
    "--supercell",
    type=int,
    default=2,
    help="provides size of supercell for calculation, default %(default)s",
)

cli.add_argument(
    "--device",
    type=str,
    default="cpu",
    choices=["cuda", "cpu", "mps"],
    help=" provide device to run. default %(default)s",
)

cli.add_argument(
    "--t_min",
    type=float,
    default=0,
    help="start temperature for CV calculations, default %(default)s",
)

cli.add_argument(
    "--t_max",
    type=float,
    default=1000,
    help="end temperature for CV calculations, default %(default)s",
)

cli.add_argument(
    "--t_step",
    type=float,
    default=50,
    help="temperature step for CV calculations, default %(default)s",
)

cli.add_argument(
    "--y_min",
    type=float,
    default=0,
    help="min y for band structure, default %(default)s",
)

cli.add_argument(
    "--y_max",
    type=float,
    default=10,
    help="max y for band structure, default %(default)s",
)

cli.add_argument(
    "--grid",
    nargs=3,
    type=int,
    default=[10, 10, 10],
    help="max y for band structure, default %(default)s",
)


args = cli.parse_args()
precision = args.precision
device = args.device
model = args.model
d3 = args.d3
dos = args.dos
fmax = args.minimize_fmax
name = args.system
ad = args.displacement
nc = args.supercell
t_min = args.t_min
t_step = args.t_step
t_max = args.t_max
ymax = args.y_max
ymin = args.y_min
m = args.grid

sysname = pathlib.PurePosixPath(name).stem

opt_log = "optimisation.log"


def optimize(a: Atoms, fmax: Float, optimize_cell: Bool, output: Str) -> None:
    x = os.getpid()
    traj = f"tmp{x}.traj"
    if optimize_cell:
        ua = FrechetCellFilter(a, hydrostatic_strain=True)
        # dyn = PreconFIRE(ua,precon=Exp(A=3), use_armijo=True,trajectory=traj,logfile=opt_log)
        dyn = FIRE(ua, trajectory=traj, logfile=opt_log)
    else:
        # dyn = PreconFIRE(a,precon=Exp(A=3), use_armijo=True,trajectory=traj,logfile=opt_log)
        dyn = FIRE(a, trajectory=traj, logfile=opt_log)

    dyn.run(fmax=fmax)
    write(output, a, format="extxyz")
    x = read(traj, index=":")
    stem = pathlib.PurePosixPath(output).stem
    write(f"{stem}-traj.xyz", x, format="extxyz")
    write(f"{stem}.cif", a, format="cif")
    try:
        p = pathlib.Path(traj)
        p.unlink()
    except:
        pass


# no magnetic moments considered
def Phonopy2ASEAtoms(s: PhonopyAtoms) -> Atoms:

    return Atoms(
        symbols=s.symbols,
        scaled_positions=s.scaled_positions,
        cell=s.cell,
        masses=s.masses,
        pbc=True,
    )


def ASE2PhonopyAtoms(s: Atoms) -> PhonopyAtoms:

    return PhonopyAtoms(
        symbols=s.get_chemical_symbols(),
        cell=s.cell.array,
        scaled_positions=s.get_scaled_positions(),
        masses=s.get_masses(),
    )


def calc_forces(calculator: Calculator, s: PhonopyAtoms) -> ArrayLike:
    atoms = Phonopy2ASEAtoms(s)
    atoms.calc = calculator
    return atoms.get_forces()


# path = [[[0, 0, 0], [0.5, 0, 0.5], [0.625, 0.25, 0.625]],
#        [[0.375, 0.375, 0.75], [0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.25, 0.75]]]
# labels = ["$\\Gamma$", "X", "U", "K", "$\\Gamma$", "L", "W"]
# qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=51)
result = {}

system = read(name)


supercell_matrix = ((nc, 0, 0), (0, nc, 0), (0, 0, nc))


mace_calc = mace_mp(device=device, model=model, default_dtype=precision, dispersion=d3)
system.calc = mace_calc
if d3:
    model = f"{model}-d3"

t = clock.datetime.now()
print(f"optimisation started at {t}")
opt_log = f"{sysname}-{model}-opt.log"
optimize(system, fmax, optimize_cell=True, output=f"{sysname}-{model}-opt.cif")
dt = clock.datetime.now() - t
print(f"phonons optimization finished after {dt}")

cell = ASE2PhonopyAtoms(system)

t = clock.datetime.now()
print(f"phonons displacements started at {t}")
phonon = phonons.Phonopy(cell, supercell_matrix)
phonon.generate_displacements(distance=ad)
disp_supercells = phonon.supercells_with_displacements
dt = clock.datetime.now() - t
print(f"phonons displacements finished after {dt}")
t = clock.datetime.now()
print(f"phonons forces started at {t}")
phonon.forces = [
    calc_forces(mace_calc, supercell)
    for supercell in disp_supercells
    if supercell is not None
]
dt = clock.datetime.now() - t
print(f"phonons forces finished after {dt}")
t = clock.datetime.now()
print(f"phonons force constants started at {t}")
phonon.produce_force_constants()
dt = clock.datetime.now() - t
print(f"phonons force constants finished after {dt}")
t = clock.datetime.now()
print(f"phonons thermal started at {t}")
phonon.run_mesh()
phonon.run_thermal_properties(t_step=t_step, t_max=t_max, t_min=t_min)
result["phonon"] = phonon
result["thermal_properties"] = phonon.get_thermal_properties_dict()
dt = clock.datetime.now() - t
print(f"phonons thermal finished after {dt}")

t = clock.datetime.now()
print(f"phonons save started at {t}")
result["phonon"].save(f"{sysname}-ase-{model}.yaml", settings={"force_constants": True})
dt = clock.datetime.now() - t
print(f"phonons save finished after {dt}")

t = clock.datetime.now()
print(f"phonons band auto started at {t}")
result["phonon"].auto_band_structure(
    write_yaml=True, filename=f"{sysname}-auto-band-{model}.yaml"
)
bplt = result["phonon"].plot_band_structure()
bplt.savefig(f"{sysname}-bs-auto-ase-{model}.pdf")

np = result["phonon"].band_structure.path_connections.count(False)
fig, ax = plt.subplots(ncols=np, nrows=1, layout="constrained")
bs = result["phonon"].band_structure.plot(ax)
_ = [ax[i].set_ylim(ymin, ymax) for i in range(np)]
_ = [ax[i].get_yaxis().set_visible(False) for i in range(1, np, 1)]
fig.savefig(f"{sysname}-bs-ase-auto-zoom-{model}.pdf")


dt = clock.datetime.now() - t
print(f"phonons band auto finished after {dt}")

t = clock.datetime.now()
print(f"phonons save cv started at {t}")
out = open(f"{sysname}-cv-{model}.dat", "w")
temp = result["thermal_properties"]["temperatures"]
cv = result["thermal_properties"]["heat_capacity"]
entropy = result["thermal_properties"]["entropy"]
free = result["thermal_properties"]["free_energy"]


print("#Temperature [K] | Cv  | H | S ", file=out)
for x, y, z, tt in zip(temp, cv, free, entropy):
    print(f"{x} {y} {z} {tt}", file=out)
out.close()
dt = clock.datetime.now() - t
print(f"phonons save cv&co finished after {dt}")

if dos:
    t = clock.datetime.now()
    print(f"phonons dos started at {t}")
    result["phonon"].run_mesh(m)
    result["phonon"].run_total_dos()
    bplt = result["phonon"].plot_total_dos()
    bplt.savefig(f"{sysname}-dos-{model}.pdf")
    dt = clock.datetime.now() - t
    print(f"phonons dos finished after {dt}")

    t = clock.datetime.now()
    print(f"phonons pdos started at {t}")
    result["phonon"].run_mesh(m, with_eigenvectors=True, is_mesh_symmetry=False)
    result["phonon"].run_projected_dos()
    bplt = result["phonon"].plot_projected_dos()
    bplt.savefig(f"{sysname}-pdos-{model}.pdf")
    dt = clock.datetime.now() - t
    print(f"phonons pdos finished after {dt}")
