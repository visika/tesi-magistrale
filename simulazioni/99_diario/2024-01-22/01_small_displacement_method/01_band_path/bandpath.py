import marimo

__generated_with = "0.1.79"
app = marimo.App()


@app.cell
def __():
    from ase.io import read
    from ase.visualize import view
    return read, view


@app.cell
def __(read):
    atoms = read("POSCAR")
    return atoms,


@app.cell
def __(atoms):
    lat = atoms.cell.get_bravais_lattice()
    return lat,


@app.cell
def __(lat):
    lat.description()
    return


@app.cell
def __(lat):
    lat.plot_bz()
    return


if __name__ == "__main__":
    app.run()
