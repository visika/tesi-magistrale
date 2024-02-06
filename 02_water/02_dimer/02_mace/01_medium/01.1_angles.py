import marimo

__generated_with = "0.1.81"
app = marimo.App()


@app.cell
def __():
    import marimo as mo
    from ase.io import read
    return mo, read


@app.cell
def __(read):
    atoms = read("01_relax-positions/final.pdb")
    return atoms,


@app.cell
def __():
    from ase.visualize import view
    return view,


@app.cell
def __(atoms, view):
    view(atoms)
    return


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        """
        Ottieni gli indici selezionando la voce di menu:
        `View > Show Labels > Atom Index`

        - 0 ossigeno accettore
        - 1 idrogeno accettore
        - 2 idrogeno accettore
        - 3 ossigeno donore
        - 4 idrogeno donore esterno (f)
        - 5 idrogeno donore che punta (d)
        """
    )
    return


@app.cell
def __(atoms):
    # Angolo alpha = 5.5°
    atoms.get_angle(0, 3, 5).round(1)
    return


@app.cell
def __(atoms):
    # Angolo interno della molecola accettore = 104.87°
    atoms.get_angle(1, 0, 2).round(2)
    return


@app.cell
def __(atoms):
    # Angolo interno della molecola donore = 104.83°
    atoms.get_angle(4, 3, 5).round(2)
    return


@app.cell
def __(atoms):
    # Distanza O-O = 291.2 pm
    atoms.get_distance(0, 3).round(2)
    return


@app.cell
def __(atoms):
    # Angolo beta = 124.4°
    atoms.get_angle(1, 0, 3).round(1)
    return


if __name__ == "__main__":
    app.run()
