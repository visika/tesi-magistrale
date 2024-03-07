#!/usr/bin/env python3

# Perform geometry optimization, keeping the cell fixed

# I/O
from ase.io import read, write

# I/O to save dynamics trajectory
from ase.io.trajectory import Trajectory

from ase.calculators.lammpslib import LAMMPSlib

import os

# Optimizer
from ase.optimize import QuasiNewton

# The CellFilter object is a wrapper that creates an Atoms object with a variable cell
from ase.constraints import UnitCellFilter

# For easier unit conversions
from ase import units

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

# Define the path to store the trajectory file
path = "relax-full/"

# Check if the directory "relas-full" exists.
# If not, create it

if not os.path.exists(str(path)):
    os.makedirs(str(path))

for pressure in [0, 1]:
    # Create new Atoms object
    geo = atoms.copy()
    # Attach calculator
    geo.calc = LAMMPScalc

    # Define the atoms object with variable cell
    ucf = UnitCellFilter(geo, scalar_pressure=pressure * units.GPa)
    # Optimizer is defined for the UnitCellFilter, not for the initial atoms object
    qn = QuasiNewton(ucf)

    # Define the Trajectory object
    traj = Trajectory(path + "simulation_" + str(pressure) + "GPa.traj", "w", geo)
    # Attach trajectory to optimizer
    qn.attach(traj)

    qn.run(fmax=0.03)
    write(
        images=geo,
        filename=path + "final_P" + str(pressure) + "GPa.pdb",
        format="proteindatabank",
    )
