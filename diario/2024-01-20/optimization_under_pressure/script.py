#!/usr/bin/env python3
import os
import numpy as np
from ase.constraints import UnitCellFilter
from ase.io import Trajectory, read, write
from ase import units
from ase.calculators.lammpslib import LAMMPSlib
from ase.optimize import QuasiNewton

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

path = "relax-full/"
os.makedirs(path, exist_ok=True)

for pressure in np.logspace(-3, 1, 10):
    geo = atoms.copy()
    geo.calc = LAMMPScalc
    ucf = UnitCellFilter(geo, scalar_pressure=pressure * units.GPa)
    qn = QuasiNewton(ucf)
    traj = Trajectory(path + "simulation_P" + f"{pressure:.2f}" + "GPa.traj", "w", geo)
    qn.attach(traj)
    qn.run(fmax=0.001)
    write(
        images=geo,
        filename=path + "final_P" + f"{pressure:.2f}" + "GPa.pdb",
        format="proteindatabank",
    )
