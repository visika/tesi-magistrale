import marimo

__generated_with = "0.1.81"
app = marimo.App()


@app.cell
def __():
    from ase.io import read
    return read,


@app.cell
def __(read):
    atoms = read("final.pdb")
    return atoms,


@app.cell
def __():
    from ase.visualize import view
    return view,


@app.cell
def __(atoms, view):
    view(atoms)
    return


@app.cell
def __():
    # 0 ossigeno accettore
    # 1 idrogeno accettore
    # 2 idrogeno accettore
    # 3 ossigeno donore
    # 5 idrogeno che punta
    return


@app.cell
def __(atoms):
    # alpha = 5.5°
    atoms.get_angle(0, 3, 5)
    return


@app.cell
def __(atoms):
    # Molecola accettore
    atoms.get_angle(1, 0, 2)
    return


@app.cell
def __(atoms):
    # Distanza O-O = 291.2 pm
    atoms.get_distance(0, 3)
    return


@app.cell
def __(atoms):
    # beta = 124.4°
    atoms.get_angle(1, 0, 3)
    return


@app.cell
def __(atoms):
    atoms.get_dihedral(0, 1, 3, 0)
    return


if __name__ == "__main__":
    app.run()
