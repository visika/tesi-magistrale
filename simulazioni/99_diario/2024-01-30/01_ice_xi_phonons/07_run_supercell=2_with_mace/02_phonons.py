#!/usr/bin/env python3
from ase.io import read
from ase.phonons import Phonons
from mace.calculators import mace_mp
import matplotlib.pyplot as plt


path = "relax-positions/"
atoms = read(path + "final.pdb")
calc = mace_mp(model="medium", dispersion=False, default_dtype="float64", device="cuda")
atoms.calc = calc

# Supercell number
N = 2
K = 10
for d in [0.05, 0.01, 1e-3, 1e-4, 1e-5]:
    print(f"[I] Phonon calculation for displacement: {d}")
    ph = Phonons(atoms, calc, supercell=(N, N, N), delta=d)
    ph.run()
    ph.read(acoustic=True)
    ph.clean()
    bandpath = atoms.cell.bandpath(npoints=100)
    bs = ph.get_band_structure(bandpath)
    dos = ph.get_dos(kpts=(K, K, K)).sample_grid(npts=100, width=1e-3)

    fig = plt.figure(1, figsize=(7, 5))
    ax = fig.add_axes([0.12, 0.07, 0.67, 0.85])
    plt.title(f"${N} \\times {N} \\times {N}$ supercell, $\\delta = {d}$, fmax=0.001")
    emin = -0.04
    emax = 0.06
    bs.plot(ax=ax, emin=emin, emax=emax)
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
    plt.savefig(f"bandstructure_d={d}.png")
    plt.close(fig)
