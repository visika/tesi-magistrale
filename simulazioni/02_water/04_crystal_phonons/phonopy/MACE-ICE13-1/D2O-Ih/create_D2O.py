import marimo

__generated_with = "0.3.3"
app = marimo.App()


@app.cell
def __():
    from ase.io import read, write
    return read, write


@app.cell
def __(read):
    atoms = read("POSCAR")
    print("Simboli:", atoms.symbols)

    print("Masse iniziali")
    print(atoms.get_masses())

    atoms.set_masses([2.01410177811 if atom.symbol == "H" else None for atom in atoms])

    print("Masse dopo ridefinizione")
    print(atoms.get_masses())
    return atoms,


@app.cell
def __(atoms, write):
    write(filename="D2O-Ih.xyz", images=atoms)
    return


@app.cell
def __(read):
    d2o = read("D2O-Ih.xyz")
    return d2o,


@app.cell
def __(d2o):
    d2o.get_masses()
    return


if __name__ == "__main__":
    app.run()
