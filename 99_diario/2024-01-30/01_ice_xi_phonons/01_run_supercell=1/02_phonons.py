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

path = "relax-positions/"
atoms = read(path + "final.pdb")
atoms.calc = LAMMPScalc

# Supercell number
N = 1
ph = Phonons(atoms, LAMMPScalc, supercell=(N, N, N), delta=0.05)
ph.run()
