import marimo

__generated_with = "0.5.0"
app = marimo.App()


@app.cell
def __():
    from ase.io import Trajectory
    import matplotlib.pyplot as plt
    import marimo as mo
    import ase.geometry.analysis
    import numpy as np
    return Trajectory, ase, mo, np, plt


@app.cell
def __(Trajectory):
    traj = Trajectory("MD-NVT/molecular_dynamics.traj", mode="r")
    return traj,


@app.cell
def __(traj):
    temp = [atoms.get_temperature() for atoms in traj]
    return temp,


@app.cell
def __(mo, plt, temp):
    plt.axhline(y=297.15, color="r", linestyle="--", label="297.15 K")
    plt.plot(temp)
    mo.mpl.interactive(plt.gcf())
    return


@app.cell
def __(traj):
    thermalized = traj[100]
    mask = [symbol == "O" for symbol in thermalized.get_chemical_symbols()]
    oxygens = thermalized[mask]
    return mask, oxygens, thermalized


@app.cell
def __(ase, mask, traj):
    rdf = []

    for atoms in traj[100:]:
        geo = atoms.copy()
        geo = geo[mask]
        data = ase.geometry.analysis.get_rdf(atoms=geo, rmax=7, nbins=40)
        rdf.append(data)
    return atoms, data, geo, rdf


@app.cell
def __(np, rdf):
    rdf_average = np.mean([e[0] for e in rdf], axis=0)
    rdf_std = np.std([e[0] for e in rdf], axis=0)
    return rdf_average, rdf_std


@app.cell
def __(plt, rdf, rdf_average):
    plt.plot(rdf[0][1], rdf_average)
    return


if __name__ == "__main__":
    app.run()
