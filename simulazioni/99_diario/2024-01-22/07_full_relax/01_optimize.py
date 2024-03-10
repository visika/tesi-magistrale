#!/usr/bin/env python3
import os
from ase.io import read, write
from ase.optimize import BFGS
from ase.calculators.lammpslib import LAMMPSlib
from ase.constraints import UnitCellFilter
from ase import units

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
pressure = 0.0
ucf = UnitCellFilter(atoms, scalar_pressure=pressure * units.GPa)
opt = BFGS(
    ucf,
    logfile=path + "optimization.log",
    trajectory=path + "optimization.traj",
)
opt.run(fmax=0.001, steps=1000)
write(images=atoms, filename=path + "final.pdb", format="proteindatabank")