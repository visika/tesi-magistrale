#!/usr/bin/env python3

# Perform geometry optimization, keeping the cell fixed

# I/O
from ase.io import read, write

# Optimizer
from ase.optimize import BFGS

# I/O to save dynamics trajectory
from ase.io.trajectory import Trajectory

# To write a good POSCAR file
from ase.build import sort

from ase.calculators.lammpslib import LAMMPSlib

import os

atoms = read("POSCAR")
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
atoms.calc = LAMMPScalc

# Define path to store the trajectory file
path = "01_relax-positions/"

# Check if directory "relax-positions" exists
# If not, create it
if not os.path.exists(str(path)):
    os.makedirs(str(path))

# Run optimization
# Create the optimizer object associated to the atoms "atoms" with calculator "LAMMPScalc"
dyn = BFGS(atoms)

# Define Trajectory object
traj = Trajectory(path + "simulation.traj", "w", atoms)

# Attach Trajectory to the optimizer
dyn.attach(traj)

# Relax the atoms positions until the forces are < 0.03 ev/Å
dyn.run(fmax=0.03)

# Write the final optimized geometry
write(images=atoms, filename=path + "final.pbd", format="proteindatabank")

write(images=sort(atoms), filename=path + "POSCAR", format="vasp")
