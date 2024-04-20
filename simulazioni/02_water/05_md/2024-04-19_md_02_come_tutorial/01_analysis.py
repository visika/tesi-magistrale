import marimo

__generated_with = "0.4.1"
app = marimo.App()


@app.cell
def __():
    import marimo as mo
    from ase import units
    from ase.io import read, Trajectory
    from ase.visualize import view
    return Trajectory, mo, read, units, view


@app.cell
def __(Trajectory):
    trajectory = Trajectory("MD-NVT/molecular_dynamics.traj", mode="r")
    return trajectory,


@app.cell
def __(trajectory, view):
    view(trajectory)
    return


@app.cell
def __(trajectory):
    image = trajectory[0]
    # Select only oxygen atoms
    mask = [s == "O" for s in image.get_chemical_symbols()]
    oxygens = image[mask]
    return image, mask, oxygens


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        """
        # Coefficiente di diffusione

        Consulta la [documentazione](https://wiki.fysik.dtu.dk/ase/ase/md.html#ase.md.analysis.DiffusionCoefficient) di ASE.
        """
    )
    return


@app.cell
def __():
    from ase.md.analysis import DiffusionCoefficient
    return DiffusionCoefficient,


@app.cell
def __(DiffusionCoefficient, trajectory, units):
    coef = DiffusionCoefficient(
        traj=trajectory,
        timestep=10 * units.fs,
    )
    return coef,


@app.cell
def __(coef):
    coef.print_data()
    return


@app.cell
def __(mo):
    mo.md(
        """
        # Radial Distribution Function

        Può essere utile usare il parametro `elements` di `get_rdf` per separare l'analisi in due.
        """
    )
    return


@app.cell
def __():
    from ase.geometry import analysis
    return analysis,


@app.cell
def __(analysis, mask, trajectory):
    analysis.get_rdf(trajectory[0][mask], rmax=7, nbins=80)
    return


@app.cell
def __(analysis, trajectory):
    analysis.get_rdf(trajectory[0], rmax=7, nbins=80, elements=[1, 8])
    return


@app.cell
def __(analysis, mask, trajectory):
    rdf = []

    rmax = 7
    nbins = 80

    for t in trajectory:
        geo = t.copy()
        # Prendi solo gli ossigeni
        geo = geo[mask]
        data = analysis.get_rdf(atoms=geo, rmax=rmax, nbins=nbins)
        rdf.append(data[0])

    # Average of the Radial Distribution Function over the images
    import numpy as np

    rdf_average = np.mean(rdf, axis=0)
    return data, geo, nbins, np, rdf, rdf_average, rmax, t


@app.cell
def __(rdf):
    rdf
    return


@app.cell
def __(rdf_average):
    rdf_average
    return


@app.cell
def __(np):
    exp = np.loadtxt("../rdf-exp", usecols=(0, 1))
    return exp,


@app.cell
def __(exp, nbins, np, rdf_average, rmax):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    x = np.arange(0, rmax, rmax / nbins)
    ax.plot(x, rdf_average, label="MACE-ICE13-1")
    plt.xlabel("r [Å]")
    plt.ylabel("$g_\mathrm{OO}(r)$")

    ax.plot(exp[:, 0], exp[:, 1], linestyle="--", label="Experiment")

    plt.grid(ls="--", alpha=0.5)
    plt.legend()
    plt.title("Radial Distribution Function (1.0 ps)")
    return ax, fig, plt, x


if __name__ == "__main__":
    app.run()
