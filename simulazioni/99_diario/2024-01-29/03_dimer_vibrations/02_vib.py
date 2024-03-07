#!/usr/bin/env python3
from ase.io import read
from ase.calculators.lammpslib import LAMMPSlib
from ase.vibrations import Vibrations

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

atoms = read("final.pdb")
atoms.center(vacuum=3.0)
atoms.calc = LAMMPScalc

for d in [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]:
    vib = Vibrations(atoms, delta=d, name=f"vib_delta={d}")
    vib.run()
    vib.summary(log=f"dimer_delta={d}_summary.txt")
