import marimo

__generated_with = "0.4.2"
app = marimo.App()


@app.cell
def __():
    import marimo as mo
    from ase import units
    return mo, units


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

    # Questa simulazione è stata fatta a 273.15 K = 0° C
    trajectory_10ps = Trajectory(
        filename="/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/05_md/2024-04-20_allunga_tempo/MD-NVT/molecular_dynamics.traj",
        mode="r",
    )

    # Simulazione NVT con la temperatura corretta di 297.15 K
    # con timestep di 0.5 fs, loginterval=10, 10000 step
    # quindi simulati in totale 5 ps
    trajectory_297K = Trajectory(
        filename="/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/05_md/2024-04-21_temperatura_corretta/MD-NVT/molecular_dynamics.traj",
        mode="r",
    )
    return trajectory_05ps, trajectory_10ps, trajectory_1ps, trajectory_297K


@app.cell
def __(Trajectory):
    # Simulazione NVT con la temperatura corretta di 297.15 K
    # con timestep di 0.5 fs, loginterval=10, 10000 step
    # quindi simulati in totale 5 ps
    trajectory_297K_MP0 = Trajectory(
        filename="/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/05_md/2024-04-22_temperatura_corretta_mace-mp-0/MD-NVT/molecular_dynamics.traj",
        mode="r",
    )
    return trajectory_297K_MP0,


@app.cell
def __(Trajectory):
    # Simulazione NVT a 297.15 K,
    # timestep di 0.5 fs, loginterval di 10, 10000 step,
    # modello MACE-MP-0 senza dispersione
    trajectory_297K_MP0_no_disp = Trajectory(
        filename="/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/05_md/2024-04-22_temperatura_corretta_mace-mp-0_senza_dispersione/MD-NVT/molecular_dynamics.traj",
        mode="r",
    )
    return trajectory_297K_MP0_no_disp,


@app.cell(disabled=True)
def __(trajectory_1ps, view):
    view(trajectory_1ps)
    return


@app.cell
def __(mo):
    mo.md("## Coefficiente di diffusione")
    return


@app.cell
def __(trajectory_1ps, units):
    from ase.md.analysis import DiffusionCoefficient

    coef = DiffusionCoefficient(
        traj=trajectory_1ps,
        timestep=10 * units.fs,
        molecule=True
    )
    return DiffusionCoefficient, coef


@app.cell
def __(coef):
    coef.print_data()
    return


@app.cell(disabled=True)
def __(DiffusionCoefficient, trajectory_1ps, units):
    coef_atoms = DiffusionCoefficient(
        traj=trajectory_1ps,
        timestep=10 * units.fs,
        molecule=False
    )

    coef_atoms.print_data()
    return coef_atoms,


@app.cell(disabled=True)
def __(DiffusionCoefficient, trajectory_10ps, units):
    # Vediamo il coefficiente di diffusione del ghiaccio a 0°C
    coef_zeroC = DiffusionCoefficient(
        traj=trajectory_10ps[200:],
        timestep=0.5*10 * units.fs,
    )

    coef_zeroC.print_data()
    return coef_zeroC,


@app.cell(disabled=True)
def __(DiffusionCoefficient, trajectory_297K, units):
    # Coefficiente di diffusione per la simulazione con MACE-ICE13-1 a 297.15 K
    # Per il timestep del DiffusionCoefficient bisogna moltiplicare il timestep della dinamica molecolare per il loginterval
    _md_timestep = 0.5 * units.fs
    _md_loginterval = 10
    coef_297K = DiffusionCoefficient(
        traj=trajectory_297K[100:],
        timestep=_md_timestep * _md_loginterval,
    )
    coef_297K.print_data()
    return coef_297K,


@app.cell
def __(DiffusionCoefficient, trajectory_297K_MP0, units):
    # Coefficiente di diffusione per la simulazione con MACE-MP-0 a 297.15 K
    _md_timestep = 0.5 * units.fs
    _md_loginterval = 10
    coef_297K_MP0 = DiffusionCoefficient(
        traj=trajectory_297K_MP0[100:],
        timestep=_md_timestep * _md_loginterval,
    )
    coef_297K_MP0.print_data()
    return coef_297K_MP0,


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
    plt.title("Temperatura del sistema,\ntempo di simulazione 1ps, termostato di Langevin 300K")
    plt.gca()
    return temperatures,


@app.cell
def __(mo):
    mo.md("Così non va bene, devo rendere più lunga la simulazione, in modo da poter scartare termini che non sono termalizzati, e trattenere abbastanza dati validi.")
    return


@app.cell
def __(plt, trajectory_10ps):
    temperatures_10ps = [a.get_temperature() for a in trajectory_10ps]
    plt.plot(temperatures_10ps, label="MD NVT")

    # Draw a horizontal line
    plt.axhline(y=273.15, color="r", linestyle="--", label="273.15 K")

    plt.grid(alpha=0.5, linestyle="--")
    plt.title("Temperatura del sistema,\ntempo di simulazione 10ps, termostato di Langevin 273.15K")
    plt.legend()
    plt.gca()
    return temperatures_10ps,


@app.cell
def __(plt, temperatures_10ps):
    # Scarto i primi 200 passi
    plt.plot(temperatures_10ps[200:], label="MD NVT")
    plt.ylim(0, 350)
    plt.grid()
    plt.gca()
    return


@app.cell
def __(plt, trajectory_297K):
    temperatures_297K = [a.get_temperature() for a in trajectory_297K]
    plt.plot(temperatures_297K, label="MD NVT")
    plt.xlabel("Step")
    plt.ylabel("Temperature [K]")

    # Draw a horizontal line
    plt.axhline(y=297.15, color="r", linestyle="--", label="297.15 K")
    plt.axvline(x=100, color="g", linestyle="-.", label="100 steps")

    plt.grid(alpha=0.5, linestyle="--")
    plt.title("Temperatura del sistema, MACE-ICE13-1\ntempo di simulazione 5 ps, termostato di Langevin 297.15K")
    plt.legend()
    plt.gca()
    return temperatures_297K,


@app.cell
def __(plt, trajectory_297K_MP0):
    temperatures_297K_MP0 = [a.get_temperature() for a in trajectory_297K_MP0]
    plt.plot(temperatures_297K_MP0, label="MD NVT")
    plt.xlabel("Step")
    plt.ylabel("Temperature [K]")

    # Draw a horizontal line
    plt.axhline(y=297.15, color="r", linestyle="--", label="297.15 K")
    plt.axvline(x=100, color="g", linestyle="-.", label="100 steps")

    plt.grid(alpha=0.5, linestyle="--")
    plt.title("Temperatura del sistema, MACE-MP-0+D3\ntempo di simulazione 5 ps, termostato di Langevin 297.15K")
    plt.legend()
    plt.gca()
    return temperatures_297K_MP0,


@app.cell
def __(plt, trajectory_297K_MP0_no_disp):
    temperatures_297K_MP0_no_disp = [a.get_temperature() for a in trajectory_297K_MP0_no_disp]
    plt.plot(temperatures_297K_MP0_no_disp, label="MD NVT")
    plt.xlabel("Step")
    plt.ylabel("Temperature [K]")

    # Draw a horizontal line
    plt.axhline(y=297.15, color="r", linestyle="--", label="297.15 K")
    plt.axvline(x=100, color="g", linestyle="-.", label="100 steps")

    plt.grid(alpha=0.5, linestyle="--")
    plt.title("Temperatura del sistema, MACE-MP-0\ntempo di simulazione 5 ps, termostato di Langevin 297.15K")
    plt.legend()
    plt.gca()
    return temperatures_297K_MP0_no_disp,


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
def __(analysis, nbins, rmax, trajectory_297K):
    _image = trajectory_297K[100]
    # Select only oxygen atoms
    _mask = [s == "O" for s in _image.get_chemical_symbols()]
    _oxygens = _image[_mask]

    rdf_297K = []

    # Considera solo le immagini dalla 100 in poi,
    # per scartare la parte di termalizzazione
    for _t in trajectory_297K[100:]:
        _geo = _t.copy()
        # Prendi solo gli ossigeni
        _geo = _geo[_mask]
        _data = analysis.get_rdf(atoms=_geo, rmax=rmax, nbins=nbins)
        rdf_297K.append(_data)
    return rdf_297K,


@app.cell
def __(np, rdf_297K):
    rdf_297K_average = np.mean([e[0] for e in rdf_297K], axis=0)
    return rdf_297K_average,


@app.cell
def __(np, rdf_297K):
    rdf_297K_std = np.std([e[0] for e in rdf_297K], axis=0)
    return rdf_297K_std,


@app.cell
def __(analysis, nbins, rmax, trajectory_297K):
    rdf_oh_297K = []
    for _t in trajectory_297K[100:]:
        _geo = _t.copy()
        _data = analysis.get_rdf(
            atoms=_geo, rmax=rmax, nbins=nbins, elements=[8, 1]
        )
        rdf_oh_297K.append(_data)
    return rdf_oh_297K,


@app.cell
def __(np, rdf_oh_297K):
    rdf_oh_297K_average = np.mean([e[0] for e in rdf_oh_297K], axis=0)
    return rdf_oh_297K_average,


@app.cell
def __(np, rdf_oh_297K):
    rdf_oh_297K_std = np.std([e[0] for e in rdf_oh_297K], axis=0)
    return rdf_oh_297K_std,


@app.cell
def __(analysis, nbins, np, rmax, trajectory_297K):
    _image = trajectory_297K[100]
    # Select only oxygen atoms
    _mask = [s == "H" for s in _image.get_chemical_symbols()]
    _hydrogens = _image[_mask]

    rdf_hh_297K = []

    # Considera solo le immagini dalla 100 in poi,
    # per scartare la parte di termalizzazione
    for _t in trajectory_297K[100:]:
        _geo = _t.copy()
        # Prendi solo gli ossigeni
        _geo = _geo[_mask]
        _data = analysis.get_rdf(atoms=_geo, rmax=rmax, nbins=nbins)
        rdf_hh_297K.append(_data)

    rdf_hh_297K_average = np.mean([e[0] for e in rdf_hh_297K], axis=0)
    rdf_hh_297K_std = np.std([e[0] for e in rdf_hh_297K], axis=0)
    return rdf_hh_297K, rdf_hh_297K_average, rdf_hh_297K_std


@app.cell
def __(analysis, nbins, rmax, trajectory_297K_MP0):
    _image = trajectory_297K_MP0[100]
    # Select only oxygen atoms
    _mask = [s == "O" for s in _image.get_chemical_symbols()]
    _oxygens = _image[_mask]

    rdf_297K_MP0 = []

    # Considera solo le immagini dalla 100 in poi,
    # per scartare la parte di termalizzazione
    for _t in trajectory_297K_MP0[100:]:
        _geo = _t.copy()
        # Prendi solo gli ossigeni
        _geo = _geo[_mask]
        _data = analysis.get_rdf(atoms=_geo, rmax=rmax, nbins=nbins)
        rdf_297K_MP0.append(_data)
    return rdf_297K_MP0,


@app.cell
def __(np, rdf_297K_MP0):
    rdf_297K_MP0_average = np.mean([e[0] for e in rdf_297K_MP0], axis=0)
    rdf_297K_MP0_std = np.std([e[0] for e in rdf_297K_MP0], axis=0)
    return rdf_297K_MP0_average, rdf_297K_MP0_std


@app.cell
def __(analysis, nbins, np, rmax, trajectory_297K_MP0_no_disp):
    _image = trajectory_297K_MP0_no_disp[100]
    # Select only oxygen atoms
    _mask = [s == "O" for s in _image.get_chemical_symbols()]
    _oxygens = _image[_mask]

    rdf_297K_MP0_no_disp = []

    # Considera solo le immagini dalla 100 in poi,
    # per scartare la parte di termalizzazione
    for _t in trajectory_297K_MP0_no_disp[100:]:
        _geo = _t.copy()
        # Prendi solo gli ossigeni
        _geo = _geo[_mask]
        _data = analysis.get_rdf(atoms=_geo, rmax=rmax, nbins=nbins)
        rdf_297K_MP0_no_disp.append(_data)

    rdf_297K_MP0_no_disp_average = np.mean(
        [e[0] for e in rdf_297K_MP0_no_disp], axis=0
    )
    rdf_297K_MP0_no_disp_std = np.std([e[0] for e in rdf_297K_MP0_no_disp], axis=0)
    return (
        rdf_297K_MP0_no_disp,
        rdf_297K_MP0_no_disp_average,
        rdf_297K_MP0_no_disp_std,
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
def __(
    nbins,
    np,
    plt,
    rdf_05ps_average,
    rdf_1ps_average,
    rdf_297K,
    rdf_297K_average,
    rmax,
):
    fig, ax = plt.subplots()

    x = np.arange(0, rmax, rmax / nbins)
    ax.plot(x, rdf_1ps_average, label="MACE-ICE13-1 1ps")
    ax.plot(x, rdf_05ps_average, label="MACE-ICE13-1 0.5ps")
    ax.plot(rdf_297K[0][1], rdf_297K_average, label="MACE-ICE13-1 297.15K")

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
def __(df, plt, rdf_297K, rdf_297K_average, rdf_297K_std, rmax):
    plt.plot(df["r"], df["g_OO(r)"], label="Experiment", color="tab:blue")
    # Plot the RDF of the simulation with standard deviation
    _asse_x = rdf_297K[0][1]
    plt.plot(
        _asse_x, rdf_297K_average, label="MACE-ICE13-1", color="tab:green"
    )
    plt.fill_between(
        _asse_x,
        rdf_297K_average - rdf_297K_std,
        rdf_297K_average + rdf_297K_std,
        alpha=0.3,
        color="tab:green",
    )

    plt.xlabel("r [Å]")
    plt.ylabel("$g_\mathrm{OO}(r)$")
    plt.grid(ls="--", alpha=0.5)

    # Set xlim to rmax
    plt.xlim(0, rmax)
    plt.legend()
    plt.title("Radial Distribution Function of liquid water\nLangevin NVT MD, T=297.15 K, simulation time: 5 ps")

    # plt.savefig("rdf.png")
    plt.gca()
    return


@app.cell
def __(
    df,
    plt,
    rdf_297K_MP0,
    rdf_297K_MP0_average,
    rdf_297K_MP0_std,
    rmax,
):
    plt.plot(df["r"], df["g_OO(r)"], label="Experiment", color="tab:blue")
    # Plot the RDF of the simulation with standard deviation
    _asse_x = rdf_297K_MP0[0][1]
    plt.plot(
        _asse_x, rdf_297K_MP0_average, label="MACE-MP-0+D3", color="tab:green"
    )
    plt.fill_between(
        _asse_x,
        rdf_297K_MP0_average - rdf_297K_MP0_std,
        rdf_297K_MP0_average + rdf_297K_MP0_std,
        alpha=0.3,
        color="tab:green",
    )

    plt.xlabel("r [Å]")
    plt.ylabel("$g_\mathrm{OO}(r)$")
    plt.grid(ls="--", alpha=0.5)

    # Set xlim to rmax
    plt.xlim(0, rmax)
    plt.legend()
    plt.title("Radial Distribution Function of liquid water\nLangevin NVT MD, T=297.15 K, simulation time: 5 ps")

    # plt.savefig("rdf.png")
    plt.gca()
    return


@app.cell
def __(
    df,
    plt,
    rdf_297K_MP0_no_disp,
    rdf_297K_MP0_no_disp_average,
    rdf_297K_MP0_no_disp_std,
    rmax,
):
    plt.plot(df["r"], df["g_OO(r)"], label="Experiment", color="tab:blue")
    # Plot the RDF of the simulation with standard deviation
    _asse_x = rdf_297K_MP0_no_disp[0][1]
    plt.plot(
        _asse_x,
        rdf_297K_MP0_no_disp_average,
        label="MACE-MP-0",
        color="tab:green",
    )
    plt.fill_between(
        _asse_x,
        rdf_297K_MP0_no_disp_average - rdf_297K_MP0_no_disp_std,
        rdf_297K_MP0_no_disp_average + rdf_297K_MP0_no_disp_std,
        alpha=0.3,
        color="tab:green",
    )

    plt.xlabel("r [Å]")
    plt.ylabel("$g_\mathrm{OO}(r)$")
    plt.grid(ls="--", alpha=0.5)

    # Set xlim to rmax
    plt.xlim(0, rmax)
    plt.legend()
    plt.title(
        "Radial Distribution Function of liquid water\nLangevin NVT MD, T=297.15 K, simulation time: 5 ps"
    )

    # plt.savefig("rdf.png")
    plt.gca()
    return


@app.cell
def __(mo):
    mo.md(
        """
        ## RDF OH

        I valori delle RDF parziali sono sottostime, perché sono scalati di un fattore.
        """
    )
    return


@app.cell
def __(plt, rdf_oh_297K, rdf_oh_297K_average, rdf_oh_297K_std):
    # plt.plot(df["r"], df["g_OO(r)"], label="Experiment", color="tab:blue")
    # Plot the RDF of the simulation with standard deviation
    _asse_x = rdf_oh_297K[0][1]
    plt.plot(_asse_x, rdf_oh_297K_average, label="MACE-ICE13-1", color="tab:green")
    plt.fill_between(
        _asse_x,
        rdf_oh_297K_average - rdf_oh_297K_std,
        rdf_oh_297K_average + rdf_oh_297K_std,
        alpha=0.3,
        color="tab:green",
    )

    plt.xlabel("r [Å]")
    plt.ylabel("$g_\mathrm{OH}(r)$")
    plt.grid(ls="--", alpha=0.5)

    # Set xlim to rmax
    plt.xlim(1.5, 5.5)
    plt.ylim(0, 2)
    plt.legend()
    plt.title(
        "Radial Distribution Function of liquid water\nLangevin NVT MD, T=297.15 K, simulation time: 5 ps"
    )

    # plt.savefig("rdf.png")
    plt.gca()
    return


@app.cell
def __(mo):
    mo.md("## RDF HH")
    return


@app.cell
def __(plt, rdf_hh_297K, rdf_hh_297K_average, rdf_hh_297K_std):
    # plt.plot(df["r"], df["g_OO(r)"], label="Experiment", color="tab:blue")
    # Plot the RDF of the simulation with standard deviation
    _asse_x = rdf_hh_297K[0][1]
    plt.plot(_asse_x, rdf_hh_297K_average, label="MACE-ICE13-1", color="tab:green")
    plt.fill_between(
        _asse_x,
        rdf_hh_297K_average - rdf_hh_297K_std,
        rdf_hh_297K_average + rdf_hh_297K_std,
        alpha=0.3,
        color="tab:green",
    )

    plt.xlabel("r [Å]")
    plt.ylabel("$g_\mathrm{HH}(r)$")
    plt.grid(ls="--", alpha=0.5)

    # Set xlim to rmax
    plt.xlim(1, 5)
    # plt.ylim(0, 2)
    plt.legend()
    plt.title(
        "Radial Distribution Function of liquid water\nLangevin NVT MD, T=297.15 K, simulation time: 5 ps"
    )

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


@app.cell
def __(mo):
    mo.md("# Densità")
    return


@app.cell(disabled=True)
def __(trajectory_297K, view):
    view([trajectory_297K[0], trajectory_297K[-1]])
    return


@app.cell
def __(trajectory_297K):
    trajectory_297K[0].get_volume()
    return


@app.cell
def __(trajectory_297K):
    trajectory_297K[-1].get_volume()
    return


@app.cell
def __(trajectory_297K):
    trajectory_297K[0].pbc
    return


@app.cell
def __(mo, trajectory_297K):
    _geo = trajectory_297K[-1].copy()
    # _geo.set_scaled_positions(trajectory_297K[-1].get_scaled_positions())
    # If coordinate is outside the unit cell, delete the atom
    _to_delete = []
    for i, pos in enumerate(_geo.get_scaled_positions(wrap=False)):
        if any(pos < 0) or any(pos > 1):
            _to_delete.append(i)
    print(f"Eliminati: {_to_delete}")
    del _geo[_to_delete]
    # view(_geo)
    _mass = _geo.get_masses().sum()
    _density = _mass / _geo.get_volume()
    mo.md(f"{_density * 1.66053906660:.2f} g/cm³")
    return i, pos


@app.cell
def __(mo):
    mo.md(
        """
        ```
        ❯ pint-convert 0.49992553545145046u/angstrom^3 g/cm^3 
        0.49992553545145046 unified_atomic_mass_unit / angstrom ** 3 = 0.83014588201(25) g/cm³
        ```
        """
    )
    return


@app.cell(disabled=True)
def __(trajectory_297K, view):
    view(trajectory_297K[-1])
    return


@app.cell
def __(mo, trajectory_297K_MP0):
    _geo = trajectory_297K_MP0[-1].copy()
    # _geo.set_scaled_positions(trajectory_297K[-1].get_scaled_positions())
    # If coordinate is outside the unit cell, delete the atom
    _to_delete = []
    for _i, _pos in enumerate(_geo.get_scaled_positions(wrap=False)):
        if any(_pos < 0) or any(_pos > 1):
            _to_delete.append(_i)
    print(f"Eliminati: {_to_delete}")
    del _geo[_to_delete]
    # view(_geo)
    _mass = _geo.get_masses().sum()
    _density = _mass / _geo.get_volume()
    mo.md(f"{_density * 1.66053906660:.2f} g/cm³")
    return


@app.cell
def __(mo, trajectory_10ps):
    # Traiettoria a 0° C
    _geo = trajectory_10ps[-1].copy()
    # _geo.set_scaled_positions(trajectory_297K[-1].get_scaled_positions())
    # If coordinate is outside the unit cell, delete the atom
    _to_delete = []
    for _i, _pos in enumerate(_geo.get_scaled_positions(wrap=False)):
        if any(_pos < 0) or any(_pos > 1):
            _to_delete.append(_i)
    print(f"Eliminati: {_to_delete}")
    del _geo[_to_delete]
    # view(_geo)
    _mass = _geo.get_masses().sum()
    _density = _mass / _geo.get_volume()
    mo.md(f"{_density * 1.66053906660:.2f} g/cm³")
    return


@app.cell
def __():
    # ❯ pint-convert 1u/angstrom^3 g/cm^3 
    # 1.0 unified_atomic_mass_unit / angstrom ** 3 = 1.66053906660(50) g/cm³
    return


if __name__ == "__main__":
    app.run()
