#!/usr/bin/env python3
# Questo script, dato un file POSCAR, ottimizza la geometria e trova le vibrazioni del cristallo.
import os
from ase.io import read, write
from ase.optimize import BFGS
from ase.phonons import Phonons
from mace.calculators import mace_mp
import matplotlib.pyplot as plt
import datetime as clock


def optimize(a, fmax: float, output_path: str) -> None:
    opt = BFGS(
        atoms, logfile=path + "optimization.log", trajectory=path + "optimization.traj"
    )
    opt.run(fmax=fmax, steps=1000)
    write(images=atoms, filename=output_path + "/final.xyz")

def plot_the_figure() -> None:
    fig = plt.figure(1, figsize=(7, 5))
    ax = fig.add_axes([0.12, 0.07, 0.67, 0.85])
    plt.title(
        f"MACE-MP-0-d3-medium: ICE Ih${N}^3$ supercell, ${K}^3$ k-points $\\delta = ${delta}, fmax={fmax}"
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
    plt.savefig(f"bandstructure_delta={delta}_N={N}.png")
    plt.close(fig)


fmax = 1e-8

atoms = read("POSCAR")
calc = mace_mp(model="medium", dispersion=True, default_dtype="float64", device="cpu")
atoms.calc = calc

path = "relax-positions/"
os.makedirs(path, exist_ok=True)

# optimize(a=atoms, fmax=fmax, output_path=path)


# Supercell number
N = 6
K = 10
delta = 1e-5

t = clock.datetime.now()
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=delta)
ph.run()
ph.read(acoustic=True)
ph.clean()
dt = clock.datetime.now() - t
print(f"phonons calculated after {dt}")

bandpath = atoms.cell.bandpath(npoints=100)
bandstructure = ph.get_band_structure(bandpath)
dos = ph.get_dos(kpts=(K, K, K)).sample_grid(npts=100, width=1e-3)

plot_the_figure()
