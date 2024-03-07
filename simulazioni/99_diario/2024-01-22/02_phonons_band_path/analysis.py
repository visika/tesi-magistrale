import marimo

__generated_with = "0.1.79"
app = marimo.App()


@app.cell
def __():
    from ase.io import read
    from ase.phonons import Phonons
    return Phonons, read


@app.cell
def __(read):
    atoms = read("POSCAR")
    return atoms,


@app.cell
def __(Phonons, atoms):
    N = 2
    ph = Phonons(atoms, supercell=(N, N, N), delta=0.05)
    ph.read(acoustic=True)
    return N, ph


@app.cell
def __(atoms):
    path = atoms.cell.bandpath(npoints=100)
    return path,


@app.cell
def __(path, ph):
    bs = ph.get_band_structure(path)
    K = 10
    dos = ph.get_dos(kpts=(K, K, K)).sample_grid(npts=100, width=1e-3)
    return K, bs, dos


@app.cell
def __(bs, dos):
    import matplotlib.pyplot as plt

    fig = plt.figure(1, figsize=(7, 4))
    ax = fig.add_axes([0.12, 0.07, 0.67, 0.85])
    emin = -0.05
    emax = 0.45
    bs.plot(ax=ax, emin=emin, emax=emax)
    dosax = fig.add_axes([0.8, 0.07, 0.17, 0.85])
    dosax.fill_between(
        dos.get_weights(),
        dos.get_energies(),
        y2=0,
        color="grey",
        edgecolor="k",
        lw=1,
    )
    dosax.set_ylim(emin, emax)
    dosax.set_yticks([])
    dosax.set_xticks([])
    dosax.set_xlabel("DOS", fontsize=16)
    plt.show()
    return ax, dosax, emax, emin, fig, plt


if __name__ == "__main__":
    app.run()
