#!/usr/bin/env python3
import os
from ase.io import write, read
from ase.calculators.lammpslib import LAMMPSlib
from ase.optimize import BFGS

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

# atoms = Atoms("HOH", positions=[(0, 0, 0), (0, 1, 0), (0, 0, 1)])
atoms = read("not_optimized.xyz")
atoms.center(vacuum=3.0)
atoms.calc = LAMMPScalc

path = "01_relax-positions/"
os.makedirs(path, exist_ok=True)

opt = BFGS(
    atoms,
    logfile=path + "optimization.log",
    trajectory=path + "optimization.traj",
    maxstep=1000,
)
opt.run(fmax=0.0001)
write(path + "final.pdb", atoms)
