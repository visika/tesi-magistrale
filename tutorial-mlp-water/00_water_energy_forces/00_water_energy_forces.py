# ase I/O
from ase.io import read, write
# LAMMPS calculator
from ase.calculators.lammpslib import LAMMPSlib

# Inizializza LAMMPS
atom_types = {'O' : 2, 'H' : 1}
cmd_nn = [      "#!/bin/bash",
                "variable runnerDir       string \"/lustre/home/tccourse/ZEN/TEST/train_008\"",
                "variable runnerCutoff    equal  8.466835984",
                "pair_style hdnnp ${runnerCutoff} dir ${runnerDir} showew no showewsum 1 resetew yes maxew 200000 cflength 1.889726 cfenergy 0.036749",
                "pair_coeff * * H O",
        ]
LAMMPScalc = LAMMPSlib(atom_types=atom_types, lmpcmds=cmd_nn, tmp_dir="./", log_file='log.lammps', keep_alive=True)

# Load a geometry and compute its energy

# Initialize the geometry
atoms = read("POSCAR")
# Attach a calculator to the atoms object
atoms.calc = LAMMPScalc

# Compute the potential energy (total electronic energy) of the loaded geometry
energy = atoms.get_potential_energy()
print("Energy:", energy)
# Compute the forces on each atom
forces = atoms.get_forces()
print('Forces:\n', forces)
