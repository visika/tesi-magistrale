import marimo

__generated_with = "0.1.88"
app = marimo.App()


@app.cell
def __():
    from ase import Atoms
    return Atoms,


@app.cell
def __(Atoms):
    a = Atoms("CuO2", positions=[
        (0., 0., 0.,),
        (4., 0., 0.),
        (0., 4., 0.)
    ], pbc=True)
    return a,


@app.cell
def __(a):
    a.center(vacuum=15.)
    return


@app.cell
def __(a):
    from ase.visualize import view
    view(a)
    return view,


@app.cell
def __():
    from ase.build import nanotube
    return nanotube,


@app.cell
def __(nanotube):
    nanotubbo = nanotube(30, 1, length=2, symbol="C", verbose=True)
    return nanotubbo,


@app.cell
def __(nanotubbo, view):
    view(nanotubbo)
    return


if __name__ == "__main__":
    app.run()
