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

opt = BFGS(atoms, logfile="optimization.log", trajectory="optimization.traj")
opt.run(fmax=0.01)
write(images=atoms, filename="optimized.pdb", format="proteindatabank")

scales = np.arange(0.8, 1.5, 0.1)

point_zero_energy = []

for s in scales:
    atoms_scaled = atoms.copy()
    atoms_scaled.calc = LAMMPScalc
    atoms_scaled.set_cell(atoms.cell * s, scale_atoms=True)
    v = atoms_scaled.cell.volume

    # opt = BFGS(
    #     atoms_scaled,
    #     logfile=f"optimization_v={v}.log",
    #     trajectory=f"optimization_v={v}.traj",
    # )
    # opt.run(fmax=0.01)
    write(
        images=atoms_scaled, filename=f"optimized_v={v}.pdb", format="proteindatabank"
    )

    potentialenergy = atoms_scaled.get_potential_energy()
    point_zero_energy.append(potentialenergy)

    N = 4
    ph = Phonons(atoms_scaled, LAMMPScalc, supercell=(N, N, N), delta=0.02 * s)
    ph.run()
    ph.read(acoustic=True)
    phonon_energies, phonon_DOS = ph.dos(kpts=(10, 10, 10), npts=1000, delta=5e-4)
    ph.clean()

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
    nmol = atoms_scaled.get_global_number_of_atoms() / 3.0
    free[:, 1] /= nmol
    np.savetxt(f"helmholtz_v={v}", free, delimiter=",")
np.savetxt("point_zero_energy.txt", point_zero_energy)
