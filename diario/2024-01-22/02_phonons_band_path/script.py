#!/usr/bin/env python3
from ase.io import read
from ase.phonons import Phonons
from ase.calculators.lammpslib import LAMMPSlib

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

atoms = read("POSCAR")
atoms.calc = LAMMPScalc

# Supercell number
N = 2
ph = Phonons(atoms, LAMMPScalc, supercell=(N, N, N), delta=0.05)
ph.run()

ph.read(acoustic=True)
# ph.clean()

path = atoms.cell.bandpath(npoints=100)
bs = ph.get_band_structure(path)

# K-points number
K = 10
dos = ph.get_dos(kpts=(K, K, K)).sample_grid(npts=100, width=1e-3)

# Plot the band structure and DOS:
import matplotlib.pyplot as plt  # noqa

fig = plt.figure(1, figsize=(7, 4))
ax = fig.add_axes([0.12, 0.07, 0.67, 0.85])

# emax = 0.035
# bs.plot(ax=ax, emin=0.0, emax=emax)
bs.write("bs.json")

# dosax = fig.add_axes([0.8, 0.07, 0.17, 0.85])
# dosax.fill_between(
#     dos.get_weights(), dos.get_energies(), y2=0, color="grey", edgecolor="k", lw=1
# )

# # dosax.set_ylim(0, emax)
# dosax.set_yticks([])
# dosax.set_xticks([])
# dosax.set_xlabel("DOS", fontsize=18)

# fig.savefig("H2O_phonon.svg")
