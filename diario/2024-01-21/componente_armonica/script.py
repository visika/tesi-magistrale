#!/usr/bin/env python3
from ase.io import read, write
from ase.calculators.lammpslib import LAMMPSlib
from ase.optimize import BFGS
from ase.phonons import Phonons
from ase.thermochemistry import CrystalThermo
import numpy as np
import pandas as pd

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

# Ottimizzazione della geometria

opt = BFGS(atoms, logfile="optimization.log", trajectory="optimization.traj")
opt.run(fmax=0.001)
write(images=atoms, filename="optimized.pdb", format="proteindatabank")

potentialenergy = atoms.get_potential_energy()

ph = Phonons(atoms, LAMMPScalc, supercell=(1, 1, 1), delta=0.02)
ph.run()
ph.read(acoustic=True)
phonon_energies, phonon_DOS = ph.dos(kpts=(10, 10, 10), npts=1000, delta=5e-4)

# Helmholtz free energy
thermo = CrystalThermo(
    phonon_DOS=phonon_DOS,
    phonon_energies=phonon_energies,
    potentialenergy=potentialenergy,
    formula_units=1,
)

free = []
for t in np.arange(10, 300, 10):
    F = thermo.get_helmholtz_energy(temperature=t)
    free.append([t, F])
free = np.array(free)
nmol = atoms.get_global_number_of_atoms() / 3.0
free[:, 1] /= nmol
np.savetxt("helmholtz", free, delimiter=",")
