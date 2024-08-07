import marimo

__generated_with = "0.7.17"
app = marimo.App(width="medium")


@app.cell
def __():
    import numpy as np
    from ase.io import read
    from mace.calculators import mace_mp
    import subprocess
    import os
    import glob
    import re
    from natsort import natsorted, ns
    return glob, mace_mp, natsorted, np, ns, os, re, read, subprocess


@app.cell
def __(mace_mp, read, subprocess):
    # Generazione della supercella
    supercell = 2

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
    return atoms, calc, inphon, supercell


@app.cell
def __(np):
    # Leggi dal file DISP gli spostamenti consigliati

    dtype = np.dtype(
        [("label", int), ("value1", float), ("value2", float), ("value3", float)]
    )

    disp = np.loadtxt("DISP", usecols=(1, 2, 3, 4), dtype=dtype)
    return disp, dtype


@app.cell
def __(atoms, calc, disp, np):
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
    return dynmat,


@app.cell
def __(disp, dynmat):
    # Calcola le forze sugli atomi per le configurazioni in cui vengono applicati gli spostamenti

    print(f"Number of displacements: {len(disp)}")

    for _j, _d in enumerate(disp):
        print(f"Loop number: {_j + 1}")
        dynmat(_j)
    return


@app.cell(disabled=True)
def __(disp, dynmat):
    import concurrent.futures
    from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor


    def progress_indicator(future):
        print(".", end="", flush=True)


    with ProcessPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(dynmat, i) for i in range(len(disp))]
        for future in futures:
            future.add_done_callback(progress_indicator)
    return (
        ProcessPoolExecutor,
        ThreadPoolExecutor,
        concurrent,
        executor,
        future,
        futures,
        progress_indicator,
    )


@app.cell
def __(disp, glob, natsorted):
    with open("FORCES", "w") as _f:
        _f.writelines(f"{str(len(disp))}\n")
        filepaths = glob.glob("DYNMAT.*")
        filepaths = natsorted(filepaths)
        for filepath in filepaths:
            with open(filepath) as _f2:
                _content = _f2.read()
                _f.writelines(_content)
    return filepath, filepaths


@app.cell
def __(subprocess):
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
    return


@app.cell
def __():
    import matplotlib.pyplot as plt
    return plt,


@app.cell
def __():
    import pandas as pd
    return pd,


@app.cell
def __(pd):
    freq1 = pd.read_csv(
        "FREQ1", sep="[ ]+", header=None, engine="python", index_col=0
    )
    freq1
    return freq1,


@app.cell
def __(pd):
    freq2 = pd.read_csv(
        "FREQ2", sep="[ ]+", header=None, engine="python", index_col=0
    )
    freq2
    return freq2,


@app.cell
def __(pd):
    freq3 = pd.read_csv(
        "FREQ3", sep="[ ]+", header=None, engine="python", index_col=0
    )
    freq3
    return freq3,


@app.cell
def __(freq1, plt):
    plt.plot(freq1)
    plt.ylim(0, 5)
    plt.savefig("freq1_s=1.svg")
    plt.gcf()
    return


@app.cell
def __(freq2, plt):
    plt.plot(freq2)
    return


@app.cell
def __(freq3, plt):
    plt.plot(freq3)
    return


@app.cell
def __(freq1, plt):
    plt.plot(freq1)
    return


if __name__ == "__main__":
    app.run()
