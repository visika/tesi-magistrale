import marimo

__generated_with = "0.1.79"
app = marimo.App()


@app.cell
def __():
    from ase.io import read
    from ase.visualize import view
    from ase.phonons import Phonons
    import matplotlib.pyplot as plt
    return Phonons, plt, read, view


@app.cell
def __(read):
    relaxed = read("relax-positions/final.pbd", format="proteindatabank")
    # view(relaxed)
    return relaxed,


@app.cell
def __(Phonons, relaxed):
    N = 2
    ph = Phonons(relaxed, supercell=(N, N, N), delta=0.05)
    ph.read(acoustic=True)
    return N, ph


@app.cell
def __(relaxed):
    path = relaxed.cell.bandpath(npoints=100)
    return path,


@app.cell
def __(path, ph):
    bs = ph.get_band_structure(path)
    K = 10
    dos = ph.get_dos(kpts=(K, K, K)).sample_grid(npts=100, width=1e-3)
    return K, bs, dos


@app.cell
def __(bs):
    bs.energies.shape
    return


@app.cell
def __(N, bs, dos, plt):
    fig = plt.figure(1, figsize=(7, 5))
    ax = fig.add_axes([0.12, 0.07, 0.67, 0.85])
    plt.title(f"${N} \\times {N} \\times {N}$ supercell, $\\delta = 0.5$, fmax=0.001")
    emin = -0.8
    emax = 0.65
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
    # plt.savefig("bandstructure.png")
    plt.show()
    return ax, dosax, emax, emin, fig


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
