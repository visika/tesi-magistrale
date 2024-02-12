#!/usr/bin/env python3
import os
from ase import units
from ase.constraints import UnitCellFilter
from ase.io import read, write
from ase.optimize import BFGS
from ase.phonons import Phonons
import matplotlib.pyplot as plt
from mace.calculators import mace_mp

atoms = read("POSCAR")
calc = mace_mp(model="medium", dispersion=False, default_dtype="float64", device="cpu")
atoms.calc = calc

directory = "relax-full/"
os.makedirs(directory, exist_ok=True)

pressure = 0.0
geo = atoms.copy()
geo.calc = calc
ucf = UnitCellFilter(geo, scalar_pressure=pressure * units.GPa)
opt = BFGS(
    ucf,
    logfile=directory + "optimization.log",
    trajectory=directory + "optimization.traj",
)
opt.run(fmax=0.0001, steps=1000)

write(images=geo, filename=directory + "final.pdb")

# Ottieni il numero di molecole di H2O nel sistema
nmol = geo.get_global_number_of_atoms() / 3.0

# Ottieni l'energia potenziale del sistema
potentialenergy = geo.get_potential_energy()

# L'energia per molecola si ottiene dividendo per il numero di molecole
e_crys = potentialenergy / nmol

f = open("e_crys.txt", "w")
f.write(str(e_crys))
f.close()

# Phonons calculation
N = 6
K = 10
d = 0.001

ph = Phonons(geo, calc, supercell=(N, N, N), delta=d)
ph.run()
ph.read(acoustic=True)
bandpath = atoms.cell.bandpath(npoints=100)
bandstructure = ph.get_band_structure(bandpath)
dos = ph.get_dos(kpts=(K, K, K)).sample_grid(npts=100, width=1e-3)

fig = plt.figure(1, figsize=(7, 5))
ax = fig.add_axes((0.12, 0.07, 0.67, 0.85))
plt.title("ICE II")
emin = -0.01
emax = 0.05
bandstructure.plot(ax=ax, emin=emin, emax=emax)

dosax = fig.add_axes((0.8, 0.07, 0.17, 0.85))
dosax.fill_between(
    dos.get_weights(),
    dos.get_energies(),
    y2=0,
    color="grey",
    edgecolor="k",
    lw=1,
)
dosax.set_ylim(emin, emax)
dosax.set_yticks([])
dosax.set_xticks([])
dosax.set_xlabel("DOS", fontsize=16)
plt.savefig(f"bandstructure_N={N}_d={d}.png")
plt.close(fig)
