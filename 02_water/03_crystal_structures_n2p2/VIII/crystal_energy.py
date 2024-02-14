#!/usr/bin/env python3
import os
from ase import units
from ase.constraints import UnitCellFilter
from ase.io import read, write
from ase.optimize import BFGS
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
calc = LAMMPSlib(
    atom_types=atom_types,
    lmpcmds=cmd_nn,
    tmp_dir="./",
    log_file="log.lammps",
    keep_alive=True,
)

atoms = read("POSCAR")
atoms.calc = calc

directory = "relax-full/"
os.makedirs(directory, exist_ok=True)

pressure = 0.0
geo = atoms.copy()
geo.calc = calc
ucf = UnitCellFilter(geo, scalar_pressure=pressure * units.GPa)
opt = BFGS(
    ucf,
    logfile=directory + "optimization.log",
    trajectory=directory + "optimization.traj",
)
opt.run(fmax=0.0001, steps=1000)

write(images=geo, filename=directory + "final.pdb")

# Ottieni il numero di molecole di H2O nel sistema
nmol = geo.get_global_number_of_atoms() / 3.0

# Ottieni l'energia potenziale del sistema
potentialenergy = geo.get_potential_energy()

# L'energia per molecola si ottiene dividendo per il numero di molecole
e_crys = potentialenergy / nmol

f = open("e_crys.txt", "w")
f.write(str(e_crys))
f.close()
