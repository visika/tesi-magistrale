#!/usr/bin/env python3
# Liquid water NVT (Nose-Hoover thermostat)

# I/O
from ase.io import read

from ase.calculators.lammpslib import LAMMPSlib

# For easier unit conversions
from ase import units

from ase.md.langevin import Langevin

import os

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

geo = read("h2o-128.pdb")
geo.calc = LAMMPScalc

# Define path to store the trajectory file
path = "MD-NVT/"

# Check if the directory "MD-NVT" exists.
# If not, create it.
if not os.path.exists(str(path)):
    os.makedirs(str(path))

dyn = Langevin(
    geo,
    timestep=0.5 * units.fs,  # MD time step in femtoseconds
    temperature_K=300.0,
    friction=0.01 / units.fs,
    trajectory=path + "simulation.traj",
    loginterval=10,
)

# Run 1000 MD steps: 1000 * 0.5 fs = 0.05 ps
dyn.run(1000)
