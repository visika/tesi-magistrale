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

print(atoms.get_potential_energy())

# Ottimizzazione della geometria

opt = BFGS(atoms, logfile="optimization.log", trajectory="optimization.traj")
opt.run(fmax=0.001)
write(images=atoms, filename="optimized.pdb", format="proteindatabank")

# Costruisci celle di volumi crescenti

atomi_scalati = []

for scale in np.arange(0.5, 2.0, 0.1):
    atoms_tmp = atoms.copy()
    atoms_tmp.calc = LAMMPScalc
    atoms_tmp.set_cell(atoms.cell * scale, scale_atoms=True)
    atomi_scalati.append(atoms_tmp)

volumi = [a.cell.volume for a in atomi_scalati]
energie_potenziali = [a.get_potential_energy() for a in atomi_scalati]
df = pd.DataFrame({"V": volumi, "U0": energie_potenziali})
df.to_csv("V_vs_U0.csv", index=False)
