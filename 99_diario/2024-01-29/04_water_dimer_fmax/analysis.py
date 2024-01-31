import marimo

__generated_with = "0.1.81"
app = marimo.App()


@app.cell
def __():
    from ase.io import read
    from ase.visualize import view
    return read, view


@app.cell
def __(read):
    atoms = read("01_relax-positions/final.pdb")
    return atoms,


@app.cell
def __(atoms, view):
    view(atoms)
    return


@app.cell
def __(atoms):
    atoms[0]
    return


@app.cell
def __(atoms):
    atoms[3]
    return


@app.cell
def __(atoms):
    atoms[5]
    return


@app.cell
def __(atoms):
    atoms.get_angle(0,3,5)
    return


@app.cell
def __(atoms):
    # Get angles
    from ase.geometry import get_angles
    angles = get_angles(atoms[0], atoms[3], atoms[5])
    print(angles)
    return angles, get_angles


if __name__ == "__main__":
    app.run()
