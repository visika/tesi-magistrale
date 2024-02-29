import marimo

__generated_with = "0.2.5"
app = marimo.App()


@app.cell
def __():
    from ase.io import read
    from ase.visualize import view
    import numpy as np
    return np, read, view


@app.cell
def __(read):
    atoms = read("init.xyz")
    return atoms,


@app.cell
def __():
    # view(atoms)
    return


@app.cell
def __():
    # Gli atomi di ossigeno sono il numero 0 e il numero 3
    return


@app.cell
def __(atoms):
    # Get position of atom 0
    print(atoms[0].position)
    return


@app.cell
def __(atoms):
    # Get position of atom 3
    print(atoms[3].position)
    return


@app.cell
def __(atoms):
    # Calculate distance between atom 3 and atom 0
    print(atoms.get_distance(0, 3))
    return


@app.cell
def __(atoms):
    # Get the vector between atom 3 and atom 0
    distance_oo = atoms.get_distance(0, 3)
    vector_oo = atoms.get_distance(0, 3, mic=True, vector=True)
    direction_oo = vector_oo / distance_oo
    distance_oo, vector_oo, direction_oo
    return direction_oo, distance_oo, vector_oo


@app.cell
def __(atoms, direction_oo, distance_oo, np, view):
    # Nella figura 5 di Sprik le posizioni sono tra circa 2.8Å e 6Å
    # Si può creare un array con le molecole distanziate progressivamente di più e ottenere la binding energy per ciascuna distanza

    atoms_array = []

    # Crea array equispaziato di displacement
    distances = np.linspace(2.8, 6, 10)
    displacements = distances - distance_oo
    print(displacements)

    for _displacement in displacements:
        atoms_distanced = atoms.copy()
        atoms_distanced[3].position += direction_oo * _displacement
        atoms_distanced[4].position += direction_oo * _displacement
        atoms_distanced[5].position += direction_oo * _displacement
        atoms_array.append(atoms_distanced)

    view(atoms_array)
    return atoms_array, atoms_distanced, displacements, distances


@app.cell
def __(atoms, view):
    # Apply distance constraints to the oxygen atoms
    from ase.constraints import FixBondLength
    constraint = FixBondLength(0, 3)
    atoms.set_constraint(constraint)
    view(atoms)
    return FixBondLength, constraint


if __name__ == "__main__":
    app.run()
