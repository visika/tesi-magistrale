
__generated_with = "0.7.19"

# %%
from datetime import datetime
print("Python started at", datetime.now())

# %%
# Dreams of multi-threading:

# import concurrent.futures
# from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor


# def progress_indicator(future):
#     print(".", end="", flush=True)


# with ProcessPoolExecutor(max_workers=8) as executor:
#     futures = [executor.submit(dynmat, i) for i in range(len(disp))]
#     for future in futures:
#         future.add_done_callback(progress_indicator)

# %%
print("Import starts at", (_t := datetime.now()))

# Numpy serve per effettuare il prodotto scalare
import numpy as np
# Ase serve per leggere il file delle posizioni degli atomi
from ase.io import read
# MACE è il calcolatore scelto
# from mace.calculators import mace_mp
from mace.calculators.mace import MACECalculator
# subprocess serve a eseguire phon
import subprocess
# glob serve a leggere l'elenco dei file DYNMAT
import glob
# natsort serve a ordinare in modo naturale i file DYNMAT
from natsort import natsorted
# copyfile serve per copiare il file POSCAR.Ih nel file POSCAR
from shutil import copyfile
# Progress bar per quando si calcolano le forze
from progress.bar import Bar
# marimo serve per i blocchi in markdown
import marimo as mo
# os per rimuovere i file
import os

print("Import finished after", datetime.now() - _t)

# %%
mo.md(
    r"""
    # Building the FORCES file
    From §3 of the PHON manual:

    > The ﬁrst line is the number of displacements, then, for each displacement, a line containing a number which indicates the position of the atom in the super-cell which has been moved, followed by the displacement (in crystal coordinates), followed by the forces on all the atoms in the super-cell (in units of eV/A and in Cartesian coordinates).
    > For central diﬀerences (LCENTRAL = .T., new default from version 1.38) the format of the file is the same, but there are twice as many displacements (for each displacement $u$ there is also $−u$).
    """
)

# %%
mo.md(
    r"""
    # Phonon dispersions
    From §3.2 of the PHON manual.
    """
)

# %%
# Generazione della supercella
supercell = 3

_inphon = f"""# Generated with Python

# number of ions types and masses
NTYPES = 2

# generate superlattice
LSUPER = .TRUE.
NDIM = {supercell} {supercell} {supercell}
"""

with open("INPHON", "w") as _f:
    _f.write(_inphon)

# Inserisci il file della cella originale
copyfile(src="POSCAR.Ih", dst="POSCAR")

print("Running PHON")
subprocess.run("./phon")

# Copia il file SPOSCAR in POSCAR, inserendo i simboli chimici al posto giusto,
# altrimenti ASE non sa di cosa si tratta

with open("SPOSCAR", "r") as _f:
    _contents = _f.readlines()

_contents.insert(5, "H O\n")

with open("POSCAR", "w") as _f:
    _f.writelines(_contents)

print("Initialize ASE atoms")
atoms = read("POSCAR")
# calc = mace_mp("medium", default_dtype="float64", device="cuda")
calc = MACECalculator(
    "/ibiscostorage/VirtualMatterLab/MACE-ICE13/MACE-ICE13-1.model",
    default_dtype="float64",
    device="cuda",
)
atoms.calc = calc

# Una volta che ASE ha ottenuto la geometria desiderata,
# ripristina il file POSCAR alla forma in cui PHON se lo aspettava
copyfile(src="SPOSCAR", dst="POSCAR")

# Leggi dal file DISP gli spostamenti consigliati

dtype = np.dtype(
    [("label", int), ("value1", float), ("value2", float), ("value3", float)]
)

disp = np.loadtxt("DISP", usecols=(1, 2, 3, 4), dtype=dtype)

# %%
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

# %%
# Calcola le forze sugli atomi per le configurazioni in cui vengono applicati gli spostamenti

print(f"Number of displacements: {len(disp)}")
print("Forces calculation starts at", (_t := datetime.now()))

bar = Bar(
    "Calculating forces", max=len(disp), check_tty=False, hide_cursor=False
)

for _j, _d in enumerate(disp):
    # print(f"Loop number: {_j + 1}")
    dynmat(_j)
    bar.next()

bar.finish()
print("Forces calculation finished after", datetime.now() - _t)

# TODO Qui si potrebbe saltare completamente il passaggio per i file DYNMAT,
# e scrivere direttamente nel file FORCES nel loop sopra

print("Building the FORCES file")

with open("FORCES", "w") as _f:
    _f.writelines(f"{str(len(disp))}\n")

    # Prendi i file DYNMAT con ordinamento naturale
    filepaths = glob.glob("DYNMAT.*")
    filepaths = natsorted(filepaths)

    for filepath in filepaths:
        with open(filepath, "r") as _f2:
            _content = _f2.read()
            _f.writelines(_content)

labels_gupta = ["G", "A", "K", "H", "M", "L", "G"]
bandpath_gupta = atoms.cell.bandpath(labels_gupta, density=0)
path_gupta = [bandpath_gupta.special_points.get(key) for key in labels_gupta]

np.savetxt("qi", path_gupta[:-1], fmt="%14.8f", newline=" \\\n")
np.savetxt("qf", path_gupta[1:], fmt="%14.8f", newline=" \\\n")

_inphon_0 = rf"""# Generated with Python

LSUPER = .F.

# number of ions types and masses
NTYPES = 2
MASS = 1.008 15.999

# q points section
LRECIP = .TRUE.
ND = {len(path_gupta) - 1}
NPOINTS = 51

QI = \
"""

_inphon_1 = rf"""
QF = \
"""

with open("INPHON", "w") as _f:
    _f.write(_inphon_0)

    with open("qi", "r") as _qi:
        _qi_contents = _qi.readlines()
    _f.writelines(_qi_contents)

    _f.write(_inphon_1)

    with open("qf", "r") as _qf:
        _qf_contents = _qf.readlines()
    _f.writelines(_qf_contents)

if os.path.exists("qi"):
    os.remove("qi")
if os.path.exists("qf"):
    os.remove("qf")

print("Phonon dispersions calculation started at", (_t := datetime.now()))
subprocess.run("./phon")
print("Phonon dispersions calculation finished after", datetime.now() - _t)