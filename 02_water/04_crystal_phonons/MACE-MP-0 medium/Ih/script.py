#!/usr/bin/env python3
# Questo script, dato un file POSCAR, ottimizza la geometria e trova le vibrazioni del cristallo.
import os
from ase.io import read, write
from ase.optimize import BFGS
from ase.phonons import Phonons
from mace.calculators import mace_mp
import matplotlib.pyplot as plt

fmax = 1e-8

atoms = read("POSCAR")
calc = mace_mp(model="medium", dispersion=True, default_dtype="float64", device="cpu")
atoms.calc = calc

path = "relax-positions/"
os.makedirs(path, exist_ok=True)

opt = BFGS(
    atoms, logfile=path + "optimization.log", trajectory=path + "optimization.traj"
)
opt.run(fmax=fmax, steps=1000)

write(images=atoms, filename=path + "final.xyz")

# Supercell number
N = 6
K = 10
delta = 1e-5

ph = Phonons(atoms, calc, supercell=(N, N, N), delta=delta)
ph.run()
ph.read(acoustic=True)
ph.clean()

bandpath = atoms.cell.bandpath(npoints=100)
bandstructure = ph.get_band_structure(bandpath)
dos = ph.get_dos(kpts=(K, K, K)).sample_grid(npts=100, width=1e-3)

# Plot the figure
fig = plt.figure(1, figsize=(7, 5))
ax = fig.add_axes([0.12, 0.07, 0.67, 0.85])
plt.title(
    f"MACE-MP-0 medium D: ICE Ih${N}^3$ supercell, ${K}^3$ k-points $\\delta = {delta}$, fmax={fmax}"
)
emin = -0.04
emax = 0.06
bandstructure.plot(ax=ax, emin=emin, emax=emax)
dosax = fig.add_axes([0.8, 0.07, 0.17, 0.85])
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
plt.savefig(f"bandstructure_delta={delta}.png")
plt.close(fig)
