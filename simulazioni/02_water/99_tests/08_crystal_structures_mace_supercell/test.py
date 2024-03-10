import marimo

__generated_with = "0.2.8"
app = marimo.App()


@app.cell
def __():
    from ase.io import read
    from ase.visualize import view
    return read, view


@app.cell
def __(read):
    atoms = read("Ih/POSCAR")
    return atoms,


@app.cell
def __(atoms):
    atoms
    return


@app.cell
def __(atoms, view):
    view(atoms)
    return


@app.cell
def __(atoms):
    supercell = atoms * 3
    supercell
    return supercell,


@app.cell
def __(supercell, view):
    view(supercell)
    return


@app.cell
def __(atoms):
    from ase.build import make_supercell
    asesupercell = make_supercell(atoms, (3,3,3))
    return asesupercell, make_supercell


@app.cell
def __(atoms):
    geo = atoms.copy() * 3
    return geo,


@app.cell
def __(geo, view):
    view(geo)
    return


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
