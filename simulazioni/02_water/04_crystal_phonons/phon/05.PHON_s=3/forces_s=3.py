print("Start Python")

import numpy as np
from ase.io import read
from mace.calculators import mace_mp
import subprocess
import os
import glob
import re
from natsort import natsorted, ns

print("Import complete")

# Generazione della supercella
supercell = 3

inphon = f"""# Generated with Python

# number of ions types and masses
NTYPES = 2

# generate superlattice
LSUPER = .TRUE.
NDIM = {supercell} {supercell} {supercell}
"""

with open("INPHON", "w") as _f:
    _f.write(inphon)

subprocess.run("./phon")

# Copia il file SPOSCAR in POSCAR

with open("SPOSCAR", "r") as _f:
    _contents = _f.readlines()

# Inserisci i simboli chimici, perché altrimenti ASE non sa di cosa si tratta
_contents.insert(5, "H O\n")

with open("POSCAR", "w") as _f:
    _f.writelines(_contents)

atoms = read("POSCAR")
calc = mace_mp("medium", default_dtype="float64")
atoms.calc = calc

# Leggi dal file DISP gli spostamenti consigliati

dtype = np.dtype(
    [("label", int), ("value1", float), ("value2", float), ("value3", float)]
)

disp = np.loadtxt("DISP", usecols=(1, 2, 3, 4), dtype=dtype)

def dynmat(i):
    displacement = disp[i]

    # Converti la tupla in array, così posso prenderne la parte di coordinate con lo slice
    d = list(displacement)

    _atoms_copy = atoms.copy()

    # L'indice fornito da PHON conta a partire da 1, Python a partire da 0
    _index = d[0] - 1

    # Lo spostamento consigliato da PHON è in coordinate cristalline
    _crystal_d = d[1:]
    # Bisogna fare il prodotto con i vettori di base della cella
    # TODO Quando la super-cella si ingrandisce, gli spostamenti devono tenerne conto
    _cartesian_d = np.dot(_crystal_d, _atoms_copy.cell)
    # Sposta gli atomi nella rappresentazione di ASE
    _atoms_copy[_index].position += _cartesian_d

    _atoms_copy.calc = calc

    _forces = _atoms_copy.get_forces()

    # Lo spostamento deve essere in coordinate cristalline
    _header = " ".join(map(str, d))
    # Le forze devono essere in coordinate cartesiane
    np.savetxt(
        f"DYNMAT.{i+1}", _forces, fmt="%14.8f", header=_header, comments=""
    )

    return 0

# Calcola le forze sugli atomi per le configurazioni in cui vengono applicati gli spostamenti

print(f"Number of displacements: {len(disp)}")

for _j, _d in enumerate(disp):
    print(f"Loop number: {_j + 1}")
    dynmat(_j)

with open("FORCES", "w") as _f:
    _f.writelines(f"{str(len(disp))}\n")
    filepaths = glob.glob("DYNMAT.*")
    filepaths = natsorted(filepaths)
    for filepath in filepaths:
        with open(filepath) as _f2:
            _content = _f2.read()
            _f.writelines(_content)

_inphon = """# Generated with Python

LSUPER = .F.

# number of ions types and masses
NTYPES = 2
MASS = 1.008 15.999

# q points section
LRECIP = .TRUE.
ND = 6
NPOINTS = 51
QI = 0.0        0.0        0.0 \\
     0.0        0.0        -0.5 \\
     0.333      0.333      0.0 \\
     0.333      0.333      -0.5 \\
     0.5        0.0        0.0 \\
     0.5        0.0        -0.5

QF = 0.0        0.0        -0.5 \\
     0.333      0.333      0.0 \\
     0.333      0.333      -0.5 \\
     0.5        0.0        0.0 \\
     0.5        0.0        -0.5 \\
     0.0        0.0        0.0
"""

with open("INPHON", "w") as _f:
    _f.write(_inphon)

subprocess.run("./phon")
