import marimo

__generated_with = "0.1.79"
app = marimo.App()


@app.cell
def __():
    import marimo as mo
    return mo,


@app.cell
def __(mo):
    mo.md(
        "# Espansione termica"
    )
    return


@app.cell
def __():
    from ase.io import read
    return read,


@app.cell
def __(read):
    atoms = read("POSCAR")
    return atoms,


@app.cell
def __(atoms):
    atoms
    return


@app.cell
def __(atoms):
    atoms.cell
    return


@app.cell
def __():
    import numpy as np
    return np,


@app.cell
def __(atoms, np):
    atomoni = []
    for v in np.arange(1., 2., 0.1):
        atoms_temp = atoms.copy()
        atoms_temp.set_cell(atoms.cell * v, scale_atoms=True)
        atomoni.append(atoms_temp)
    return atomoni, atoms_temp, v


@app.cell
def __(atomoni):
    atomoni
    return


@app.cell
def __(atomoni):
    [a.cell.volume for a in atomoni]
    return


@app.cell
def __(atomoni):
    from ase.visualize import view
    view(atomoni)
    return view,


@app.cell
def __(glob, re):
    import pandas as pd
    import os
    all_files = glob.glob("Helmholtz_vs_T_volume=*")
    # Sort all_files naturally
    all_files = sorted(all_files, key=lambda x: float(re.findall("volume=(.*)\.csv", x)[0]))
    all_files
    return all_files, os, pd


@app.cell
def __(all_files, pd, re):
    li = []
    for filename in all_files:
        df = pd.read_csv(filename, header=None, names=["T", "F", "V"])
        # Make volume as index
        df["V"] = float(re.findall("volume=(.*)\.csv", filename)[0])
        # Round float
        df["V"] = df["V"].round(2)
        # df.set_index("volume", inplace=True)
        # print(df)
        # Add volume as column
        li.append(df)
    li[0]
    return df, filename, li


@app.cell
def __(all_files, pd):
    df3 = pd.read_csv(all_files[0], header=None, names=["T", "F"])
    return df3,


@app.cell
def __(df3):
    # Plot the values in df3
    df3.plot(x="T", y="F")
    return


@app.cell
def __(np, pd, volumes):
    # Read values of punto_zero
    zero = pd.read_csv("punto_zero.csv", header=None, names=["U"])
    # Associate each volume to the value of punto_zero
    zero["V"] = np.round(volumes, 2)
    zero
    return zero,


@app.cell
def __(li, pd, volumes, zero):
    # Combine dataframes in li, using V as a multi-index
    df2 = pd.concat(li, axis=0, keys=volumes, names=["V", "index"])
    # df2.drop(columns=["V"], inplace=True)
    t10 = df2[df2["T"] == 10.0].copy()
    t10["U"] = zero["U"].values
    t10["tot"] = t10["U"] + t10["F"]
    t10.plot(x="V", y="tot", marker="o")
    return df2, t10


@app.cell
def __(t10):
    t10.plot(x="V", y="U", marker="o")
    return


@app.cell
def __():
    import glob
    files = glob.glob("Helmholtz_vs_T_volume=*")
    files
    return files, glob


@app.cell
def __(files):
    # Read files starting with "Helmholtz_vs_T_volume=" and extract the number after the = sign and before the file extension
    import re
    volumes = [float(re.findall("volume=(.*)\.csv", f)[0]) for f in files]
    # Sort the volumes
    volumes = sorted(volumes)
    volumes
    return re, volumes


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
