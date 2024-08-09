#!/usr/bin/env python3
#
# Questo script serve a calcolare le forze dovute agli spostamenti finiti.
# Successivamente a questo viene l'analisi delle componenti armoniche.

# nc è il lato della supercella da usare
nc = 2

import phonopy
from ase import Atoms
from datetime import datetime
from ase.optimize import BFGS
from ase.io import read, write
from mace.calculators.mace import MACECalculator
from phonopy.structure.atoms import PhonopyAtoms
from ase.md.npt import NPT
from ase.filters import FrechetCellFilter
from ase import units


def optimize(a: Atoms, fmax: float, output: str = "optimized.xyz") -> None:
    # Il filtro serve a impostare la pressione desiderata
    filter = FrechetCellFilter(atoms, scalar_pressure=1.6e3 * units.bar)
    opt = BFGS(filter, trajectory="optimization.traj", logfile="optimization.log")
    opt.run(fmax=fmax)
    write(filename=output, images=a)


def ASE2PhonopyAtoms(s: Atoms) -> PhonopyAtoms:
    return PhonopyAtoms(
        symbols=s.get_chemical_symbols(),
        cell=s.cell.array,
        scaled_positions=s.get_scaled_positions(),
        masses=s.get_masses(),
    )


def calc_forces(calculator, s: PhonopyAtoms):
    atoms = Phonopy2ASEAtoms(s)
    atoms.calc = calculator
    return atoms.get_forces()


def Phonopy2ASEAtoms(s: PhonopyAtoms) -> Atoms:
    return Atoms(
        symbols=s.symbols,
        scaled_positions=s.scaled_positions,
        cell=s.cell,
        masses=s.masses,
        pbc=True,
    )


atoms = read(filename="POSCAR")

calculator = MACECalculator(
    "/ibiscostorage/VirtualMatterLab/MACE-ICE13/MACE-ICE13-1.model",
    default_dtype="float64",
    device="cuda",
)

atoms.calc = calculator

t = datetime.now()
print(f"Optimization started at {t}")
optimize(atoms, fmax=1e-4)
dt = datetime.now() - t
print(f"Optimization finished after {dt}")

cell = ASE2PhonopyAtoms(atoms)

t = datetime.now()
print(f"Displacements started at {t}")
supercell_matrix = ((nc, 0, 0), (0, nc, 0), (0, 0, nc))
ph = phonopy.Phonopy(cell, supercell_matrix)
ph.generate_displacements(distance=0.01)
disp_supercells = ph.supercells_with_displacements
dt = datetime.now() - t
print(f"Displacements finished after {dt}")

t = datetime.now()
print(f"Forces calculation started at {t}")
ph.forces = [
    calc_forces(calculator, supercell)
    for supercell in disp_supercells
    if supercell is not None
]
dt = datetime.now() - t
print(f"Forces calculation finished after {dt}")

t = datetime.now()
print(f"Force constants calculation started at {t}")
ph.produce_force_constants()
dt = datetime.now() - t
print(f"Force constants calculation finished after {dt}")

ph.save("phonopy_params.yaml", settings={"force_constants": True})
print("Goodbye from Python!")
