import marimo

__generated_with = "0.4.2"
app = marimo.App()


@app.cell
def __():
    import marimo as mo
    return mo,


@app.cell
def __(mo):
    mo.md("# Visualizzazione dei dati sperimentali")
    return


@app.cell
def __():
    import pandas as pd
    return pd,


@app.cell
def __(pd):
    # Read a txt file, skip lines from 1 to 10 and the 12th and use tab as a separator
    # df = pd.read_csv('Ambient_water_xray_data.txt', sep='\t')
    df = pd.read_csv(
        "Ambient_water_xray_data.txt",
        sep="\t",
        skiprows=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11],
    )
    return df,


@app.cell
def __(df):
    df
    return


@app.cell
def __(df):
    # Plot g_OO(r) vs r
    import matplotlib.pyplot as plt
    plt.plot(df["r"], df["g_OO(r)"])
    # Add grid, slightly transparent and dashed
    plt.grid(alpha=0.5, linestyle="--")
    plt.gca()
    return plt,


@app.cell
def __(mo):
    mo.md("# Visualizzazione dei dati delle simulazioni")
    return


@app.cell
def __():
    from ase.geometry import analysis
    from ase.io import Trajectory
    return Trajectory, analysis


@app.cell
def __(Trajectory):
    trajectory_1ps = Trajectory(
        filename="/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/05_md/2024-04-19_md_02_come_tutorial/MD-NVT/molecular_dynamics.traj",
        mode="r",
    )

    trajectory_05ps = Trajectory(
        filename="/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/05_md/2024-04-19_md_01/MD-NVT/molecular_dynamics.traj",
        mode="r",
    )
    return trajectory_05ps, trajectory_1ps


@app.cell
def __(mo):
    mo.md(
        """
        ## Grafico della temperatura nel corso della simulazione

        Mi serve capire quali stati non ancora termalizzati mi conviene scartare, perché mi sporcano le medie della RDF.
        """
    )
    return


@app.cell
def __(plt, trajectory_1ps):
    temperatures = [a.get_temperature() for a in trajectory_1ps]
    plt.plot(temperatures)
    plt.grid(alpha=0.5, linestyle="--")
    plt.gca()
    return temperatures,


@app.cell
def __(mo):
    mo.md("Così non va bene, devo rendere più lunga la simulazione, in modo da poter scartare termini che non sono termalizzati, e trattenere abbastanza dati validi.")
    return


@app.cell
def __(analysis, trajectory_1ps):
    import numpy as np

    image = trajectory_1ps[0]
    # Select only oxygen atoms
    mask = [s == "O" for s in image.get_chemical_symbols()]
    oxygens = image[mask]

    rdf_1ps = []

    rmax = 7
    nbins = 80

    for t in trajectory_1ps:
        geo = t.copy()
        # Prendi solo gli ossigeni
        geo = geo[mask]
        data = analysis.get_rdf(atoms=geo, rmax=rmax, nbins=nbins)
        rdf_1ps.append(data[0])

    # Average of the Radial Distribution Function over the images
    rdf_1ps_average = np.mean(rdf_1ps, axis=0)
    # Standard deviation of the Radial Distribution Function over the images
    rdf_1ps_std = np.std(rdf_1ps, axis=0)
    return (
        data,
        geo,
        image,
        mask,
        nbins,
        np,
        oxygens,
        rdf_1ps,
        rdf_1ps_average,
        rdf_1ps_std,
        rmax,
        t,
    )


@app.cell
def __(analysis, nbins, np, rmax, trajectory_05ps):
    _image = trajectory_05ps[0]
    # Select only oxygen atoms
    _mask = [s == "O" for s in _image.get_chemical_symbols()]
    _oxygens = _image[_mask]

    rdf_05ps = []

    for _t in trajectory_05ps:
        _geo = _t.copy()
        # Prendi solo gli ossigeni
        _geo = _geo[_mask]
        _data = analysis.get_rdf(atoms=_geo, rmax=rmax, nbins=nbins)
        rdf_05ps.append(_data[0])

    # Average of the Radial Distribution Function over the images
    rdf_05ps_average = np.mean(rdf_05ps, axis=0)
    # Standard deviation of the Radial Distribution Function over the images
    rdf_05ps_std = np.std(rdf_05ps, axis=0)
    return rdf_05ps, rdf_05ps_average, rdf_05ps_std


@app.cell
def __(nbins, np, plt, rdf_05ps_average, rdf_1ps_average, rmax):
    fig, ax = plt.subplots()

    x = np.arange(0, rmax, rmax / nbins)
    ax.plot(x, rdf_1ps_average, label="MACE-ICE13-1 1ps")
    ax.plot(x, rdf_05ps_average, label="MACE-ICE13-1 0.5ps")

    # ax.plot(exp[:, 0], exp[:, 1], linestyle="--", label="Experiment")

    plt.xlabel("r [Å]")
    plt.ylabel("$g_\mathrm{OO}(r)$")
    plt.grid(ls="--", alpha=0.5)
    plt.legend()
    plt.title("Radial Distribution Function")
    return ax, fig, x


@app.cell
def __(mo):
    mo.md("## RDF della prima e dell'ultima istanza")
    return


@app.cell
def __(analysis, nbins, np, plt, rmax, trajectory_1ps):
    _image = trajectory_1ps[0]
    _image_last = trajectory_1ps[-1]
    # Select only oxygen atoms
    _mask = [s == "O" for s in _image.get_chemical_symbols()]
    _oxygens = _image[_mask]
    _oxygens_last = _image_last[_mask]

    rdf = analysis.get_rdf(atoms=_oxygens, rmax=rmax, nbins=nbins)[0]
    rdf_last = analysis.get_rdf(atoms=_oxygens_last, rmax=rmax, nbins=nbins)[0]

    _fig, _ax = plt.subplots()

    _x = np.arange(0, rmax, rmax / nbins)
    _ax.plot(_x, rdf, label="MACE-ICE13-1 0ps")
    _ax.plot(_x, rdf_last, label="MACE-ICE13-1 1ps")

    plt.legend()
    return rdf, rdf_last


@app.cell
def __(mo):
    mo.md("# Confronto di simulazione e dati sperimentali")
    return


@app.cell
def __(df, plt, rdf_1ps_average, rdf_1ps_std, rmax, x):
    plt.plot(df["r"], df["g_OO(r)"], label="Experiment", color="tab:blue")
    # Plot the RDF of the simulation with standard deviation
    plt.plot(x, rdf_1ps_average, label="MACE-ICE13-1 1ps", color="tab:green")
    plt.fill_between(
        x,
        rdf_1ps_average - rdf_1ps_std,
        rdf_1ps_average + rdf_1ps_std,
        alpha=0.3,
        color="tab:green",
    )

    plt.xlabel("r [Å]")
    plt.ylabel("$g_\mathrm{OO}(r)$")
    plt.grid(ls="--", alpha=0.5)

    # Set xlim to rmax
    plt.xlim(0, rmax)
    plt.legend()
    plt.title("Radial Distribution Function")

    # plt.savefig("rdf.png")
    plt.gca()
    return


@app.cell
def __(mo):
    mo.md("# Studio della RDF del cristallo di ghiaccio Ih")
    return


@app.cell
def __():
    from ase.io import read
    return read,


@app.cell
def __(read):
    ice = read("/home/mariano/Progetti/tesi-magistrale/strutture/ICE13/Ih/POSCAR")
    return ice,


@app.cell
def __(ice):
    # Mi conviene avere un sistema con lo stesso numero di atomi di quello che sto studiando, per poter fare un confronto diretto, cioè 128.

    # Conta il numero di atomi
    len(ice)
    return


@app.cell
def __(ice):
    # Di quanto devo moltiplicare?
    128 / len(ice)
    return


@app.cell
def __(ice):
    # Moltiplica la cella in tutte le direzioni di due
    double_ice = ice * 2
    len(double_ice)
    return double_ice,


@app.cell
def __(ice):
    triple_ice = ice * 3
    len(triple_ice)
    return triple_ice,


@app.cell
def __():
    from ase.visualize import view
    return view,


@app.cell
def __():
    # view([ice, double_ice])
    return


@app.cell
def __(analysis, double_ice, nbins, plt):
    _mask = [s == "O" for s in double_ice.get_chemical_symbols()]
    _oxygens = double_ice[_mask]

    _rdf_ice = analysis.get_rdf(atoms=_oxygens, rmax=5, nbins=nbins)

    plt.plot(_rdf_ice[1], _rdf_ice[0])
    return


@app.cell
def __(analysis, nbins, plt, rmax, triple_ice):
    _mask = [s == "O" for s in triple_ice.get_chemical_symbols()]
    _oxygens = triple_ice[_mask]

    rdf_ice = analysis.get_rdf(
        _oxygens,
        rmax=rmax,
        nbins=nbins,
    )

    plt.plot(rdf_ice[1], rdf_ice[0])
    return rdf_ice,


@app.cell
def __(df, plt, rdf_1ps_average, rdf_1ps_std, rdf_ice, rmax, x):
    _fig, _ax1 = plt.subplots(layout="constrained")

    _ax1.plot(df["r"], df["g_OO(r)"], label="Experiment", color="tab:blue")
    # Plot the RDF of the simulation with standard deviation
    _ax1.plot(x, rdf_1ps_average, label="MACE-ICE13-1 1ps", color="tab:green")
    _ax1.fill_between(
        x,
        rdf_1ps_average - rdf_1ps_std,
        rdf_1ps_average + rdf_1ps_std,
        alpha=0.3,
        color="tab:green",
    )
    _ax1.set_xlabel("r [Å]")
    _ax1.set_ylabel("$g_\mathrm{OO}(r)$ liquid water")
    _ax1.grid(ls="--", alpha=0.5)

    # Use second y axis for the RDF of the ice
    _ax2 = _ax1.twinx()
    _ax2.plot(
        rdf_ice[1], rdf_ice[0], label="Ice Ih", color="tab:red", linestyle="--"
    )
    _ax2.set_ylabel("$g_\mathrm{OO}(r)$ ice Ih")

    # Set xlim to rmax
    _ax1.set_xlim(0, rmax)
    _fig.legend()
    plt.title("Radial Distribution Function")

    # plt.savefig("rdf.png")
    plt.gca()
    return


if __name__ == "__main__":
    app.run()
