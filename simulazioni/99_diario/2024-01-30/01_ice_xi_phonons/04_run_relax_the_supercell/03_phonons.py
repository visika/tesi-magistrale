#!/usr/bin/env python3
from ase.io import read
from ase.phonons import Phonons
from ase.calculators.lammpslib import LAMMPSlib
import matplotlib.pyplot as plt

# Inizializza LAMMPS
atom_types = {"O": 2, "H": 1}
cmd_nn = [
    "#!/bin/bash",
    'variable runnerDir       string "/ibiscostorage/VirtualMatterLab/train_008"',
    "variable runnerCutoff    equal  8.466835984",
    "pair_style hdnnp ${runnerCutoff} dir ${runnerDir} showew no showewsum 1 resetew yes maxew 200000 cflength 1.889726 cfenergy 0.036749",
    "pair_coeff * * H O",
]
LAMMPScalc = LAMMPSlib(
    atom_types=atom_types,
    lmpcmds=cmd_nn,
    tmp_dir="./",
    log_file="log.lammps",
    keep_alive=True,
)

path = "relax-positions/"
atoms = read(path + "final.pdb")
atoms.calc = LAMMPScalc

# Supercell number
N = 1
K = 10
for d in [0.05, 0.01, 1e-3, 1e-4, 1e-5]:
    print(f"[I] Phonon calculation for displacement: {d}")
    ph = Phonons(atoms, LAMMPScalc, supercell=(N, N, N), delta=d)
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
