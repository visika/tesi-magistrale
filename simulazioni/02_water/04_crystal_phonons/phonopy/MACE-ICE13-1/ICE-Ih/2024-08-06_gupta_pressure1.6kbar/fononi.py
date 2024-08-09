import marimo

__generated_with = "0.7.17"
app = marimo.App(width="medium")


@app.cell
def __():
    from ase.phonons import Phonons
    from ase.io import read
    from mace.calculators import mace_mp
    import numpy as np
    from ase.thermochemistry import CrystalThermo
    import matplotlib.pyplot as plt
    return CrystalThermo, Phonons, mace_mp, np, plt, read


@app.cell
def __(mace_mp, read):
    atoms = read("POSCAR")
    # atoms.set_cell(atoms.cell * 2, scale_atoms=True)
    print(atoms.get_volume())
    calculator = mace_mp("medium")
    atoms.calc = calculator
    potentialenergy = atoms.get_potential_energy()
    print(potentialenergy)
    return atoms, calculator, potentialenergy


@app.cell
def __(Phonons, atoms, calculator):
    ph = Phonons(atoms, calculator, supercell=(1, 1, 1), delta=0.02)
    ph.run()
    ph.read(acoustic=True)
    return ph,


@app.cell
def __(ph):
    phonon_energies, phonon_DOS = ph.dos(kpts=(10, 10, 10), npts=1000, delta=5e-4)
    return phonon_DOS, phonon_energies


@app.cell
def __(CrystalThermo, phonon_DOS, phonon_energies, potentialenergy):
    thermo = CrystalThermo(
        phonon_energies=phonon_energies,
        phonon_DOS=phonon_DOS,
        potentialenergy=potentialenergy,
        formula_units=1,
    )
    return thermo,


@app.cell
def __(atoms, np, thermo):
    free = []  # where we save the free energy
    for temp in np.arange(10, 300, 10):
        F = thermo.get_helmholtz_energy(
            temperature=temp
        )  # exact free energy at a fixed temperature
        free.append([temp, F])
    free = np.array(free)

    nmol = atoms.get_global_number_of_atoms() / 3.0  # number of water molecules
    free[:, 1] /= nmol  # free energy per water molecule
    return F, free, nmol, temp


@app.cell
def __(free, np, plt):
    blue = "#009fff"

    fig, ax = plt.subplots(figsize=(7, 4))

    ax.plot(
        free[:, 0],
        free[:, 1],
        linestyle="--",
        linewidth=2.0,
        marker="o",
        mec="black",
        color=blue,
        markersize=10,
    )

    ax.tick_params(direction="out", width=1.5, length=6)

    plt.xticks(np.arange(0, 325, 25), fontsize=14)
    plt.yticks(fontsize=14)

    plt.ylabel("Free energy [eV/mol]", fontsize=14, labelpad=8)
    plt.xlabel("Temperature [K]", fontsize=14, labelpad=8)

    plt.title("Helmholtz Quasi-Harmonic Free Energy", fontsize=14)
    plt.grid(ls="--", alpha=0.3)
    plt.tight_layout()
    plt.show()
    return ax, blue, fig


@app.cell
def __(atoms, ph, plt):
    path = atoms.cell.bandpath('GA', npoints=10)
    bs = ph.get_band_structure(path)

    _fig = plt.figure(1, figsize=(7, 4))
    _ax = _fig.add_axes([0.12, 0.07, 0.67, 0.85])

    emax = 0.02
    bs.plot(ax=_ax, emin=0.0, emax=emax)
    return bs, emax, path


if __name__ == "__main__":
    app.run()
