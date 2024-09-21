import marimo

__generated_with = "0.8.18"
app = marimo.App(app_title="Analisi RDF")


@app.cell
def __():
    import marimo as mo
    from ase import units
    import matplotlib.pyplot as plt
    from ase.geometry import analysis
    from ase.io import Trajectory
    from ase.io import write
    from ase.visualize import view
    import pandas as pd
    return Trajectory, analysis, mo, pd, plt, units, view, write


@app.cell
def __(mo):
    mo.md(
        """
        # Visualizzazione dei dati sperimentali

        I dati descrivono la radial distribution function tra gli atomi di ossigeno.
        Nell'ultima colonna si trovano anche gli errori sui dati, che è possibile graficare.
        """
    )
    return


@app.cell
def __(pd):
    # Read a txt file, skip lines from 1 to 10 and the 12th and use tab as a separator
    # df = pd.read_csv('Ambient_water_xray_data.txt', sep='\t')
    df = pd.read_csv(
        "Ambient_water_xray_data.txt",
        sep="\t",
        skiprows=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11],
    )
    df
    return (df,)


@app.cell
def __(df, mo, plt):
    # Plot g_OO(r) vs r

    plt.plot(df["r"], df["g_OO(r)"])
    plt.fill_between(
        df["r"],
        df["g_OO(r)"] - df["error"],
        df["g_OO(r)"] + df["error"],
        alpha=0.3,
    )
    # Add grid, slightly transparent and dashed
    plt.grid(alpha=0.5, linestyle="--")

    plt.title("Experimental Radial Distribution Function")
    mo.mpl.interactive(plt.gcf())

    # DONE metti errore
    return


@app.cell
def __(mo):
    mo.md("""# Lettura delle traiettorie e visualizzazione dei dati delle simulazioni""")
    return


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
def __(trajectory_297K, view):
    view(trajectory_297K)
    return


@app.cell
def __(trajectory_297K, write):
    write(filename="traiettoria.pdb", images=trajectory_297K)
    return


@app.cell
def __(Trajectory):
    # Attenzione, il timestep di questa simulazione è differente, timestep=100
    traj_100ps = Trajectory(
        "2024-05-08_NVT_MACE-ICE13-1_100ps/MD-NVT/molecular_dynamics.traj"
    )
    return (traj_100ps,)


@app.cell
def __(traj_100ps):
    len(traj_100ps)
    return


@app.cell
def __(traj_100ps, view):
    view(traj_100ps[0])
    return


@app.cell
def __(traj_100ps):
    len(traj_100ps[0])/3
    return


@app.cell
def __(traj_100ps, view):
    view(traj_100ps[2000])
    return


@app.cell
def __(Trajectory):
    # Simulazione NVT con la temperatura corretta di 297.15 K
    # con timestep di 0.5 fs, loginterval=10, 10000 step
    # quindi simulati in totale 5 ps
    trajectory_297K_MP0 = Trajectory(
        filename="/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/05_md/2024-04-22_temperatura_corretta_mace-mp-0/MD-NVT/molecular_dynamics.traj",
        mode="r",
    )
    return (trajectory_297K_MP0,)


@app.cell
def __(Trajectory):
    # Simulazione NVT a 297.15 K,
    # timestep di 0.5 fs, loginterval di 10, 10000 step,
    # modello MACE-MP-0 senza dispersione
    trajectory_297K_MP0_no_disp = Trajectory(
        filename="/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/05_md/2024-04-22_temperatura_corretta_mace-mp-0_senza_dispersione/MD-NVT/molecular_dynamics.traj",
        mode="r",
    )
    return (trajectory_297K_MP0_no_disp,)


@app.cell(disabled=True)
def __(trajectory_1ps, view):
    view(trajectory_1ps)
    return


@app.cell
def __(mo):
    mo.md("""## Coefficiente di diffusione""")
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


@app.cell(disabled=True)
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
    return (coef_atoms,)


@app.cell(disabled=True)
def __(DiffusionCoefficient, trajectory_10ps, units):
    # Vediamo il coefficiente di diffusione del ghiaccio a 0°C
    coef_zeroC = DiffusionCoefficient(
        traj=trajectory_10ps[200:],
        timestep=0.5*10 * units.fs,
    )

    coef_zeroC.print_data()
    return (coef_zeroC,)


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
    return (coef_297K,)


@app.cell(disabled=True)
def __(DiffusionCoefficient, trajectory_297K_MP0, units):
    # Coefficiente di diffusione per la simulazione con MACE-MP-0 a 297.15 K
    _md_timestep = 0.5 * units.fs
    _md_loginterval = 10
    coef_297K_MP0 = DiffusionCoefficient(
        traj=trajectory_297K_MP0[100:],
        timestep=_md_timestep * _md_loginterval,
    )
    coef_297K_MP0.print_data()
    return (coef_297K_MP0,)


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
    return (temperatures,)


@app.cell
def __(mo):
    mo.md("""Così non va bene, devo rendere più lunga la simulazione, in modo da poter scartare termini che non sono termalizzati, e trattenere abbastanza dati validi.""")
    return


@app.cell
def __(mo, plt, trajectory_10ps):
    temperatures_10ps = [a.get_temperature() for a in trajectory_10ps]
    plt.plot(temperatures_10ps, label="MD NVT", linewidth=0.5)

    # Draw a horizontal line
    plt.axhline(y=273.15, color="r", linestyle="--", label="273.15 K")

    plt.grid(alpha=0.5, linestyle="--")
    plt.title("System temperature,\nsim. time 10ps, Langevin therm. 273.15K")
    plt.legend()

    plt.savefig("Grafici/temperature_NVT.png")

    mo.mpl.interactive(plt.gca())

    # TODO Effettuare simulazioni a tempi più lunghi, per l'acqua parliamo di 10-100 ps
    return (temperatures_10ps,)


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
    return (temperatures_297K,)


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
    return (temperatures_297K_MP0,)


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
    return (temperatures_297K_MP0_no_disp,)


@app.cell
def __(mo, plt):
    def grafico_temperatura(traiettoria, title=None, ymax=350):
        plt.plot([a.get_temperature() for a in traiettoria])

        plt.xlabel("Step")
        plt.ylabel("Temperature [K]")

        # Set ylim
        plt.ylim(0, ymax)
        ax = plt.gca()
        ax.spines[["top", "left", "right"]].set_visible(False)

        plt.title(title)

        return mo.mpl.interactive(plt.gcf())
    return (grafico_temperatura,)


@app.cell
def __(grafico_temperatura, traj_100ps):
    grafico_temperatura(traj_100ps, title="MACE-ICE13-1 100 ps")
    return


@app.cell
def __(mo):
    mo.md("""## Calcolo delle RDF""")
    return


@app.cell
def __(analysis):
    def calcola_rdf(traiettoria, rmax=7, nbins=40):
        image = traiettoria[100]
        mask = [symbol == "O" for symbol in image.get_chemical_symbols()]
        oxygens = image[mask]

        rdf = []

        for atoms in traiettoria[100:]:
            geo = atoms.copy()
            geo = geo[mask]
            data = analysis.get_rdf(atoms=geo, rmax=rmax, nbins=nbins)
            rdf.append(data)

        return rdf
    return (calcola_rdf,)


@app.cell
def __(np):
    def calcola_rdf_average(rdf):
        return np.mean([e[0] for e in rdf], axis=0)

    def calcola_rdf_std(rdf):
        return np.std([e[0] for e in rdf], axis=0)
    return calcola_rdf_average, calcola_rdf_std


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
    # TODO riduci nbins a 40 per verificare se si riduce la banda di deviazione std
    # TODO approfondisci autocorrelazione

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


@app.cell(disabled=True)
def __(analysis, nbins, rmax, trajectory_297K):
    rdf_oh_297K = []
    for _t in trajectory_297K[100:]:
        _geo = _t.copy()
        _data = analysis.get_rdf(
            atoms=_geo, rmax=rmax, nbins=nbins, elements=[8, 1]
        )
        rdf_oh_297K.append(_data)
    return (rdf_oh_297K,)


@app.cell
def __(np, rdf_oh_297K):
    rdf_oh_297K_average = np.mean([e[0] for e in rdf_oh_297K], axis=0)
    return (rdf_oh_297K_average,)


@app.cell
def __(np, rdf_oh_297K):
    rdf_oh_297K_std = np.std([e[0] for e in rdf_oh_297K], axis=0)
    return (rdf_oh_297K_std,)


@app.cell(disabled=True)
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
    return (rdf_297K_MP0,)


@app.cell
def __(calcola_rdf, traj_100ps):
    rdf_100ps = calcola_rdf(traj_100ps)
    return (rdf_100ps,)


@app.cell
def __(calcola_rdf_average, calcola_rdf_std, rdf_100ps):
    rdf_100ps_average = calcola_rdf_average(rdf_100ps)
    rdf_100ps_std = calcola_rdf_std(rdf_100ps)
    return rdf_100ps_average, rdf_100ps_std


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


@app.cell(hide_code=True)
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


@app.cell(hide_code=True)
def __(mo):
    mo.md("""## RDF della prima e dell'ultima istanza""")
    return


@app.cell(hide_code=True)
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


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        """
        # Confronto di simulazione e dati sperimentali
        ## RDF OO
        """
    )
    return


@app.cell(hide_code=True)
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
    mo.md("""### Traiettoria a 297.15 K con MACE-ICE13-1, tempo di simulazione 5 ps""")
    return


@app.cell
def __(mo):
    mo.md("""#### Analisi con nbins=80""")
    return


@app.cell
def __(analysis, nbins, np, rmax, trajectory_297K):
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

    rdf_297K_average = np.mean([e[0] for e in rdf_297K], axis=0)
    rdf_297K_std = np.std([e[0] for e in rdf_297K], axis=0)
    return rdf_297K, rdf_297K_average, rdf_297K_std


@app.cell
def __(df, mo, plt, rdf_297K, rdf_297K_average, rdf_297K_std, rmax):
    plt.plot(df["r"], df["g_OO(r)"], label="Experiment", color="tab:blue")

    plt.fill_between(
        df["r"],
        df["g_OO(r)"] - df["error"],
        df["g_OO(r)"] + df["error"],
        alpha=0.3,
    )

    # Plot the RDF of the simulation with standard deviation
    _asse_x = rdf_297K[0][1]
    plt.plot(_asse_x, rdf_297K_average, label="MACE-ICE13-1", color="tab:green")
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
    plt.title(
        "Radial Distribution Function of liquid water\nLangevin NVT MD, T=297.15 K, simulation time: 5 ps, nbins=80"
    )

    plt.ylim(top=3.0, bottom=-0.5)

    # plt.savefig("rdf.png")
    mo.mpl.interactive(plt.gcf())

    # TODO per calcolare D, simulazione NVE

    # TODO per calcolare densità, simulazione NPT
    return


@app.cell
def __(mo):
    mo.md("""#### Analisi con nbins=40""")
    return


@app.cell
def __(analysis, df, mo, np, plt, rmax, trajectory_297K):
    _image = trajectory_297K[100]
    # Select only oxygen atoms
    _mask = [s == "O" for s in _image.get_chemical_symbols()]
    _oxygens = _image[_mask]

    _rdf = []
    _nbins = 40

    # Considera solo le immagini dalla 100 in poi,
    # per scartare la parte di termalizzazione
    for _t in trajectory_297K[100:]:
        _geo = _t.copy()
        # Prendi solo gli ossigeni
        _geo = _geo[_mask]
        _data = analysis.get_rdf(atoms=_geo, rmax=rmax, nbins=_nbins)
        _rdf.append(_data)

    _rdf_average = np.mean([e[0] for e in _rdf], axis=0)
    _rdf_std = np.std([e[0] for e in _rdf], axis=0)

    plt.plot(df["r"], df["g_OO(r)"], label="Experiment", color="tab:blue")

    plt.fill_between(
        df["r"],
        df["g_OO(r)"] - df["error"],
        df["g_OO(r)"] + df["error"],
        alpha=0.3,
    )

    # Plot the RDF of the simulation with standard deviation
    _asse_x = _rdf[0][1]
    plt.plot(_asse_x, _rdf_average, label="MACE-ICE13-1", color="tab:green")
    plt.fill_between(
        _asse_x,
        _rdf_average - _rdf_std,
        _rdf_average + _rdf_std,
        alpha=0.3,
        color="tab:green",
    )

    plt.xlabel("r [Å]")
    plt.ylabel("$g_\mathrm{OO}(r)$")
    plt.grid(ls="--", alpha=0.5)

    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["top"].set_visible(False)

    # Set xlim to rmax
    plt.xlim(0, rmax)
    plt.ylim(top=3.0, bottom=-0.5)
    plt.legend()
    plt.title(
        f"Radial Distribution Function of liquid water\nLangevin NVT MD, T=297.15 K, simulation time: 5 ps, nbins={_nbins}"
    )

    # plt.savefig("rdf.png")
    mo.mpl.interactive(plt.gcf())
    return


@app.cell
def __(df, mo, plt, rmax):
    def grafica_rdf(
        rdf,
        rdf_average,
        rdf_std,
        model=None,
        simulation_time=None,
        temperature=None,
    ):
        plt.plot(df["r"], df["g_OO(r)"], label="Experiment", color="tab:blue")

        plt.fill_between(
            df["r"],
            df["g_OO(r)"] - df["error"],
            df["g_OO(r)"] + df["error"],
            alpha=0.3,
            color="tab:blue",
        )

        plt.plot(rdf[0][1], rdf_average, color="tab:green", label=model)
        plt.fill_between(
            rdf[0][1],
            rdf_average - rdf_std,
            rdf_average + rdf_std,
            alpha=0.3,
            color="tab:green",
        )

        plt.xlabel("r [Å]")
        plt.ylabel("$g_\mathrm{OO}(r)$")
        plt.grid(ls="--", alpha=0.5)

        plt.xlim(0, rmax)
        plt.ylim(top=3.0, bottom=-0.5)

        plt.legend()

        ax = plt.gca()
        ax.spines[["top", "bottom", "right"]].set_visible(False)

        plt.title(
            f"Radial Distribution Function of liquid water, Langevin NVT\n{model}, T={temperature} K, simulation time: {simulation_time}"
        )

        return mo.mpl.interactive(plt.gcf())
    return (grafica_rdf,)


@app.cell
def __(grafica_rdf, rdf_100ps, rdf_100ps_average, rdf_100ps_std):
    grafica_rdf(
        rdf_100ps,
        rdf_100ps_average,
        rdf_100ps_std,
        model="MACE-ICE13-1",
        simulation_time="100 ps",
        temperature=297.15,
    )
    return


@app.cell(hide_code=True)
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
    plt.title(
        "Radial Distribution Function of liquid water\nLangevin NVT MD, T=297.15 K, simulation time: 5 ps"
    )

    # plt.savefig("rdf.png")
    plt.gca()
    return


@app.cell
def __(
    df,
    mo,
    plt,
    rdf_297K_MP0_no_disp,
    rdf_297K_MP0_no_disp_average,
    rdf_297K_MP0_no_disp_std,
    rmax,
):
    plt.plot(df["r"], df["g_OO(r)"], label="Experiment", color="tab:blue")
    plt.fill_between(
        df["r"],
        df["g_OO(r)"] - df["error"],
        df["g_OO(r)"] + df["error"],
        alpha=0.3,
    )
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

    # Remove the right frame
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["top"].set_visible(False)

    # Set xlim to rmax
    plt.xlim(0, rmax)
    plt.legend()
    plt.title(
        "Radial Distribution Function of liquid water\nLangevin NVT MD, T=297.15 K, simulation time: 5 ps"
    )

    plt.savefig("Grafici/rdf_oo_mace-mp-0_NVT_T=297.15_t=5ps.svg")
    mo.mpl.interactive(plt.gca())
    return


@app.cell(hide_code=True)
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
    mo.md("""## RDF HH""")
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
    mo.md("""# Studio della RDF del cristallo di ghiaccio Ih""")
    return


@app.cell
def __():
    from ase.io import read
    return (read,)


@app.cell
def __(read):
    ice = read("/home/mariano/Progetti/tesi-magistrale/strutture/ICE13/Ih/POSCAR")
    return (ice,)


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
    return (double_ice,)


@app.cell
def __(ice):
    triple_ice = ice * 3
    len(triple_ice)
    return (triple_ice,)


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
    return (rdf_ice,)


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
    mo.md("""# Densità""")
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
def __(trajectory_297K):
    _v = trajectory_297K[-1].get_volume()
    _m = trajectory_297K[-1].get_masses().sum()
    _m / _v * 1.66053906660
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


@app.cell
def __():
    # Creo una cella con meno molecole d'acqua
    128 / 4
    return


@app.cell
def __():
    from ase import Atoms
    from ase.lattice.hexagonal import Hexagonal
    small_cell = Hexagonal(symbol="H2O")
    return Atoms, Hexagonal, small_cell


@app.cell
def __(small_cell, view):
    view(small_cell)
    return


@app.cell
def __(view):
    from ase.collections import g2
    from ase.build import molecule
    view(molecule("H2O"))
    return g2, molecule


@app.cell
def __():
    from ase.build import bulk
    return (bulk,)


@app.cell
def __():
    return


@app.cell
def __():
    from ase.lattice.hexagonal import HexagonalFactory


    class HexagonalH2OFactory(HexagonalFactory):
        bravais_basis = [
            [7.816, 11.083, 5.235],
            [8.147, 10.193, 5.372],
            [8.561, 11.602, 5.396],
        ]
        element_basis = (1, 0, 0)


    HEX_H2O = HexagonalH2OFactory()
    return HEX_H2O, HexagonalFactory, HexagonalH2OFactory


@app.cell
def __(HEX_H2O):
    lattice = HEX_H2O(
        symbol=("H", "O"), latticeconstant={"a": 1, "b": 1, "c": 1}, size=(1, 1, 1)
    )
    return (lattice,)


@app.cell
def __(lattice, view):
    view(lattice)
    return


@app.cell
def __(read):
    molecules = read(
        "/home/mariano/Progetti/tesi-magistrale/strutture/128_molecules/h2o-128.pdb"
    )
    return (molecules,)


@app.cell
def __(molecules, view):
    view(molecules)
    return


@app.cell
def __(read):
    cif = read("/home/mariano/Scaricati/H2O-Ice-Ih.cif")
    return (cif,)


@app.cell
def __(cif, view):
    view(cif)
    return


@app.cell
def __(read):
    ice_ih = read(
        "/home/mariano/Progetti/tesi-magistrale/strutture/ICE13/Ih/POSCAR"
    )
    return (ice_ih,)


@app.cell
def __(ice_ih):
    len(ice_ih) / 3
    return


@app.cell
def __(ice_ih, view):
    view(ice_ih * (2, 1, 1))
    return


@app.cell
def __(ice_ih):
    _v = ice_ih.get_volume()
    _m = ice_ih.get_masses().sum()
    _m / _v * 1.66053906660
    return


@app.cell
def __(calculate_density, molecule):
    the_molecule = molecule("H2O")
    cell_size = 3.11
    the_molecule.set_cell([cell_size] * 3)
    the_molecule.set_pbc(True)
    the_molecule.center()
    calculate_density(the_molecule)
    return cell_size, the_molecule


@app.cell
def __(the_molecule):
    the_molecule.get_cell().get_bravais_lattice()
    return


@app.cell
def __(the_molecule):
    the_molecule.pbc
    return


@app.cell
def __(the_molecule):
    the_molecules = the_molecule * 2
    return (the_molecules,)


@app.cell
def __(the_molecules, view):
    view(the_molecules)
    return


@app.cell
def __(the_molecules):
    the_molecules.get_cell().get_bravais_lattice()
    return


@app.cell
def __(calculate_density, the_molecules):
    calculate_density(the_molecules)
    return


@app.cell
def __(the_molecules):
    len(the_molecules) / 3
    return


@app.cell
def __():
    def calculate_density(atoms):
        v = atoms.get_volume()
        m = atoms.get_masses().sum()
        return m / v * 1.66053906660
    return (calculate_density,)


@app.cell
def __(mo):
    mo.md("""# Analisi degli effetti di size""")
    return


@app.cell
def __(mo):
    mo.md("""## 8 molecole""")
    return


@app.cell
def __(Trajectory):
    traj_8molecules = Trajectory(
        "2024-05-09_effetti_di_size/8_molecules/MD-NVT/molecular_dynamics.traj"
    )
    return (traj_8molecules,)


@app.cell
def __(traj_8molecules, view):
    view(traj_8molecules)
    return


@app.cell
def __(grafico_temperatura, traj_8molecules):
    grafico_temperatura(traj_8molecules, title="8 molecole", ymax=None)
    return


@app.cell
def __(calcola_rdf, traj_8molecules):
    rdf_8molecules = calcola_rdf(traj_8molecules, rmax=3)
    return (rdf_8molecules,)


@app.cell
def __(calcola_rdf_average, calcola_rdf_std, rdf_8molecules):
    rdf_8molecules_average = calcola_rdf_average(rdf_8molecules)
    rdf_8molecules_std = calcola_rdf_std(rdf_8molecules)
    return rdf_8molecules_average, rdf_8molecules_std


@app.cell
def __(
    grafica_rdf,
    rdf_8molecules,
    rdf_8molecules_average,
    rdf_8molecules_std,
):
    # Grafico 8 molecole
    grafica_rdf(
        rdf_8molecules,
        rdf_8molecules_average,
        rdf_8molecules_std,
        model="MACE-ICE13-1",
        simulation_time="10 ps",
        temperature=297.15,
    )
    return


@app.cell
def __(mo):
    mo.md("""## 27 molecole""")
    return


@app.cell
def __(Trajectory):
    traj_27molecules = Trajectory(
        "2024-05-09_effetti_di_size/27_molecules/MD-NVT/molecular_dynamics.traj"
    )
    return (traj_27molecules,)


@app.cell
def __(grafico_temperatura, traj_27molecules):
    grafico_temperatura(traj_27molecules, title="27 molecole", ymax=None)
    return


@app.cell
def __(
    calcola_rdf,
    calcola_rdf_average,
    calcola_rdf_std,
    traj_27molecules,
):
    rdf_27molecules = calcola_rdf(traj_27molecules, rmax=4.6, nbins=40)
    rdf_27molecules_average = calcola_rdf_average(rdf_27molecules)
    rdf_27molecules_std = calcola_rdf_std(rdf_27molecules)
    return rdf_27molecules, rdf_27molecules_average, rdf_27molecules_std


@app.cell
def __(
    grafica_rdf,
    rdf_27molecules,
    rdf_27molecules_average,
    rdf_27molecules_std,
):
    grafica_rdf(
        rdf_27molecules,
        rdf_27molecules_average,
        rdf_27molecules_std,
        model="MACE-ICE13-1",
        simulation_time="10 ps",
        temperature=297.15,
    )
    return


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
