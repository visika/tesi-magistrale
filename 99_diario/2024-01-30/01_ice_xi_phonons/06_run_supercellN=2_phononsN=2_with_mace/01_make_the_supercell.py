import marimo

__generated_with = "0.1.81"
app = marimo.App()


@app.cell
def __():
    from ase.io import read, write
    from ase.visualize import view
    return read, view, write


@app.cell
def __(read):
    atoms = read("POSCAR")
    return atoms,


@app.cell
def __():
    # view(atoms)
    return


@app.cell
def __(atoms):
    atoms.cell
    return


@app.cell
def __(atoms):
    atoms.get_cell()
    return


@app.cell
def __(atoms):
    atoms.set_pbc = (True, True, False)
    return


@app.cell
def __(atoms):
    supercell = atoms.repeat((2, 2, 2))
    return supercell,


@app.cell
def __(supercell, write):
    write("ice-xi_supercell=2.pdb", supercell)
    return


if __name__ == "__main__":
    app.run()
