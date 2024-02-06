#!/usr/bin/env python3
import os
from ase.io import read, write
from ase.calculators.lammpslib import LAMMPSlib
from ase.optimize import BFGS
from ase.vibrations import Vibrations

cutoffs = [8.466835984]

# Effettua i conti con diversi cutoff, per comprendere se questo parametro
# influenza le qualità di stabilità nel calcolo delle vibrazioni del dimero.

for cutoff in cutoffs:
    print("Cutoff:", cutoff)
    directory = f"cutoff={cutoff}"
    os.makedirs(directory, exist_ok=True)
    os.chdir(directory)
    # Inizializza LAMMPS
    atom_types = {"O": 2, "H": 1}
    cmd_nn = [
        "#!/bin/bash",
        'variable runnerDir       string "/ibiscostorage/VirtualMatterLab/train_008"',
        f"variable runnerCutoff    equal  {cutoff}",
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

    atoms = read("../init.xyz")
    atoms.center(vacuum=3.0)
    atoms.calc = LAMMPScalc

    opt = BFGS(atoms, logfile="optimization.log", trajectory="optimization.traj")
    opt.run(fmax=0.0001, steps=1000)
    write("final.pdb", atoms)

    alpha = atoms.get_angle(0, 3, 5).round(1)
    print("alpha:", alpha)

    beta = atoms.get_angle(1, 0, 3).round(1)
    print("beta:", beta)

    oo_distance = atoms.get_distance(0, 3).round(2)
    print("distanza OO:", oo_distance)

    angolo_h2o_a = atoms.get_angle(1, 0, 2).round(2)
    print("angolo interno molecola accettore:", angolo_h2o_a)

    angolo_h2o_d = atoms.get_angle(4, 3, 5).round(2)
    print("angolo interno molecola donore:", angolo_h2o_d)

    for d in [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]:
        vib = Vibrations(atoms, delta=d, name=f"vib_delta={d}")
        vib.run()
        vib.summary(log=f"dimer_delta={d}_summary.txt")
        vib.clean()

    os.chdir("..")
