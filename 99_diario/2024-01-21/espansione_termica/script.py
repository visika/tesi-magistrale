#!/usr/bin/env python3
from ase.io import read, write
from ase.calculators.lammpslib import LAMMPSlib
from ase.optimize import BFGS
from ase.phonons import Phonons
from ase.thermochemistry import CrystalThermo
import numpy as np

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

for scale in np.arange(1.0, 2.0, 0.1):
    atoms_temp = atoms.copy()
    atoms_temp.calc = LAMMPScalc
    atoms_temp.set_cell(atoms.cell * scale, scale_atoms=True)
    atomi_scalati.append(atoms_temp)

volumi = [a.cell.volume for a in atomi_scalati]
energie_potenziali = [a.get_potential_energy() for a in atomi_scalati]
np.savetxt("punto_zero.csv", energie_potenziali, delimiter=",")

# Calcolo dell'energia di Helmholtz con i fononi

for i, v in enumerate(volumi):
    a = atomi_scalati[i]
    potentialenergy = a.get_potential_energy()
    ph = Phonons(a, LAMMPScalc, supercell=(1, 1, 1), delta=0.02)
    ph.run()
    ph.read(acoustic=True)
    phonon_energies, phonon_DOS = ph.dos(kpts=(10, 10, 10), npts=1000, delta=5e-4)
    thermo = CrystalThermo(
        phonon_DOS=phonon_DOS,
        phonon_energies=phonon_energies,
        formula_units=1,
        potentialenergy=potentialenergy,
    )

    # Calcola l'energia di Helmholtz
    free = []
    for temperature in np.arange(10, 300, 10):
        F = thermo.get_helmholtz_energy(temperature=temperature)
        free.append([temperature, F])
    free = np.array(free)

    # Dividi per il numero di molecole, per ottenere l'energia di Helmholtz per molecola
    num_molecules = a.get_global_number_of_atoms() / 3.0
    free[:, 1] /= num_molecules

    np.savetxt(f"Helmholtz_vs_T_volume={v}.csv", free, delimiter=",")
