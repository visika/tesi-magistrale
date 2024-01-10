#!/usr/bin/env python3

# Phonons calculation

import numpy as np
import matplotlib.pyplot as plt

# I/O
from ase.io import read

from ase.calculators.lammpslib import LAMMPSlib

from ase.phonons import Phonons
from ase.thermochemistry import CrystalThermo

import os, shutil

atom_types = {"O": 2, "H": 1}
cmd_nn = [
    "#!/bin/bash",
    'variable runnerDir       string "/lustre/home/tccourse/ZEN/TEST/train_008"',
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

path = "./phonon"
if os.path.exists(path):
    shutil.rmtree(path)

# Start from a relaxed geometry
geo = read("../02_water_full_geometry_optimization/relax-full/final_P0GPa.pdb")
geo.calc = LAMMPScalc
potentialenergy = geo.get_potential_energy()

# Compute phonons with small displacement method
# Set super-cell and Δx (in Å)
ph = Phonons(geo, LAMMPScalc, supercell=(1, 1, 1), delta=0.02)
ph.run()
ph.read(acoustic=True)
phonon_energies, phonon_DOS = ph.dos(kpts=(10, 10, 10), npts=1000, delta=5e-4)

# Calculate the Helmholtz energy
thermo = CrystalThermo(
    phonon_energies=phonon_energies,
    phonon_DOS=phonon_DOS,
    potentialenergy=potentialenergy,
    formula_units=1,
)

# Use the phonons to calculate the Helmholtz energy

# Define an array into which save the Helmholtz energy
helmholtz = []

# Get exact Helmholtz energy at fixed temperatures
for temp in np.arange(10, 300, 10):
    F = thermo.get_helmholtz_energy(temperature=temp)
    helmholtz.append([temp, F])

helmholtz = np.array(helmholtz)

# Number of water molecules
nmol = geo.get_global_number_of_atoms() / 3.0

# Helmholtz energy per water molecule
helmholtz[:, 1] /= nmol

# Plot the results

blue = "#009fff"

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(
    helmholtz[:, 0],
    helmholtz[:, 1],
    linestyle="--",
    linewidth=2.0,
    marker="o",
    mec="black",
    color=blue,
    markersize=10,
)
ax.tick_params(direction="out", width=1.5, length=6)
plt.xticks(np.arange(0, 325, 25), fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel("Helmholtz energy [eV/mol]", fontname="arial", fontsize=14, labelpad=8)
plt.xlabel("Temperature [K]", fontname="arial", fontsize=14, labelpad=8)
plt.title("Helmholtz quasi-harmonic energy", fontname="arial", fontsize=14)
plt.grid(ls="--", alpha=0.3)
plt.tight_layout()
plt.savefig("plot.svg")
