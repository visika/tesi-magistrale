#!/usr/bin/env python3
from ase.io import Trajectory, write
from ase import Atoms
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

atoms = Atoms("HOH", positions=[(0, 0, 0), (0, 1, 0), (0, 0, 1)])
atoms.center(vacuum=3.0)
atoms.calc = LAMMPScalc

opt = BFGS(atoms)
traj = Trajectory("optimization.traj", "w", atoms)
opt.attach(traj)
opt.run(fmax=0.0001)
write("final.pdb", atoms)

potentialenergy = atoms.get_potential_energy()

f = open("e_gas.txt", "w")
f.write(str(potentialenergy))
f.close()
