import marimo

__generated_with = "0.3.3"
app = marimo.App()


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        """
        # Studio dei modi vibrazionali del polimorfo Ih del ghiaccio

        Si sono calcolate le forze e la matrice delle costanti di forza col metodo degli spostamenti finiti, usando una supercella \(2 \\times 2 \\times 2\).
        Queste quantità fisiche sono conservate nel file `phonopy_params.yaml`.

        Con gli strumenti di ASE è possibile visualizzare la prima zona di Brillouin della geometria in esame, ottenendo le coordinate dei punti speciali.
        """
    )
    return


@app.cell
def __():
    import phonopy
    import marimo as mo
    from ase.io import read
    from ase.visualize import view
    import matplotlib.pyplot as plt
    from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
    from phonopy.phonon.band_structure import get_band_qpoints
    return (
        get_band_qpoints,
        get_band_qpoints_and_path_connections,
        mo,
        phonopy,
        plt,
        read,
        view,
    )


@app.cell
def __(read):
    atoms = read(
        "/home/mariano/Progetti/tesi-magistrale/strutture/ICE13/Ih/POSCAR"
    )
    # view(atoms)
    bandpath = atoms.cell.bandpath()
    return atoms, bandpath


@app.cell
def __(atoms):
    atoms.cell.reciprocal()
    return


@app.cell
def __(bandpath, mo):
    mo.md(
        f"""
        Bandpath:  
        {bandpath}

        Special points:
        """
    )
    return


@app.cell
def __(bandpath):
    bandpath.plot()
    return


@app.cell
def __(mo):
    mo.md("Bisogna costruire un percorso. Scegliamo per ora di seguire quello studiato nell'articolo di Strässle, Saitta, Klotz del 2004, che segue il percorso \(A\Gamma K M \Gamma\).")
    return


@app.cell
def __(bandpath, get_band_qpoints_and_path_connections):
    labels = ["A", "G", "K", "M", "G"]
    path = [bandpath.special_points.get(key) for key in labels]

    qpoints, connections = get_band_qpoints_and_path_connections(
        # qpoints = get_band_qpoints(
        [path],
        npoints=101,
        # rec_lattice=atoms.cell.reciprocal()
    )
    return connections, labels, path, qpoints


@app.cell
def __(connections, labels, phonopy, qpoints):
    ph = phonopy.load("phonopy_params.yaml")
    ph.run_band_structure(
        qpoints, path_connections=connections, labels=labels, is_legacy_plot=False
    )
    return ph,


@app.cell
def __(mo):
    mo.md("## Plot di tutte le bande trovate")
    return


@app.cell
def __(ph):
    ph.plot_band_structure().gca()
    return


@app.cell
def __(mo):
    mo.md("## Zoom sulle bande più basse")
    return


@app.cell
def __(ph, plt):
    _fig, _ax = plt.subplots(
        ncols=2, nrows=1, layout="constrained", width_ratios=[9999, 1]
    )
    _band_structure_zoom = ph.band_structure.plot(_ax)

    _ax[0].set_ylim(0, 5)
    _ax[0].grid(axis="x")

    _fig.delaxes(_ax[1])
    _fig
    return


@app.cell
def __(mo):
    mesh = mo.ui.slider(1, 16, value=4, show_value=True).form()
    return mesh,


@app.cell
def __(mesh, mo):
    mo.vstack(
        [
            mo.md("## Calcolo della DOS"),
            mo.hstack(
                [
                    mo.md(
                        f"""
        Scegli di seguito il lato della griglia per il calcolo della DOS:
        """
                    ),
                    mesh,
                ]
            ),
        ]
    )
    return


@app.cell
def __(mesh, ph, plt):
    ph.run_mesh([mesh.value, mesh.value, mesh.value])
    ph.run_total_dos(
        sigma=None,
        freq_min=0,
        freq_max=14,
        freq_pitch=None,
        use_tetrahedron_method=True,
    )
    # ph.write_total_dos(filename=f"total_dos_mesh_{mesh.value}.dat")
    total_dos = ph.get_total_dos_dict()

    fig, ax = plt.subplots(layout="constrained")

    ax.plot(
        total_dos["frequency_points"],
        total_dos["total_dos"],
        color="red",
    )
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)
    ax.set_xlabel("Freq (THz)")
    ax.set_ylabel("DOS (1/THz)")
    ax.set_title(f"{mesh.value}x{mesh.value}x{mesh.value} mesh")


    def THz2K(THz):
        return THz * 47.9924307337


    def K2THz(K):
        return K * 0.0208366191233


    secax = ax.secondary_xaxis("top", functions=(THz2K, K2THz))
    secax.set_xlabel("E (K)")
    ax
    return K2THz, THz2K, ax, fig, secax, total_dos


@app.cell(disabled=True)
def __(mesh, ph, plt):
    ph.run_mesh([mesh.value, mesh.value, mesh.value])
    ph.run_total_dos(freq_min=0, freq_max=14)
    # ph.write_total_dos(filename=f"total_dos_mesh_{mesh.value}.dat")

    bplt = ph.plot_total_dos()
    bplt.xlim((0, 14))
    # bplt.ylim((0, 8))
    _ax = plt.gca()
    _ax
    return bplt,


@app.cell
def __(mo):
    mo.md("# Calcolo della capacità termica \(C_V\)")
    return


@app.cell
def __(mo):
    t_min = mo.ui.slider(
        start=0,
        stop=999,
        debounce=False,
        label="t_min",
        value=0,
        show_value=True,
        orientation="vertical",
    ).form()
    t_max = mo.ui.slider(
        start=1,
        stop=1000,
        debounce=False,
        label="t_max",
        value=100,
        show_value=True,
        orientation="vertical",
    ).form()
    t_step = mo.ui.slider(
        start=1,
        stop=100,
        debounce=False,
        label="t_step",
        value=10,
        show_value=True,
        orientation="vertical",
    ).form()
    return t_max, t_min, t_step


@app.cell
def __(mo, t_max, t_min, t_step):
    mo.hstack([t_min, t_max, t_step])
    return


@app.cell
def __(numero_molecole, ph, plt, t_max, t_min, t_step):
    ph.run_mesh(mesh=100.0)
    # ph.run_mesh(mesh=[_m] * 3)
    ph.run_thermal_properties(
        t_step=t_step.value, t_min=t_min.value, t_max=t_max.value
    )
    thermal_properties = ph.get_thermal_properties_dict()
    thermal_properties["heat_capacity"] = (
        thermal_properties["heat_capacity"] / numero_molecole
    )

    R_KJmol = 8.31446261815

    plt.plot(
        thermal_properties["temperatures"],
        thermal_properties["heat_capacity"] / R_KJmol,
        marker="o",
    )

    plt.xlabel("Temperature (K)")
    # Unità di misura definite in
    # https://phonopy.github.io/phonopy/setting-tags.html#thermal-properties-related-tags
    # you have to divide the value by number of formula unit in your unit cell by yourself
    plt.ylabel("Heat capacity (J/K/mol)")
    return R_KJmol, thermal_properties


@app.cell
def __(atoms):
    numero_molecole = len(atoms) / 3
    numero_molecole
    return numero_molecole,


@app.cell
def __(mo):
    mo.md(
        """
        ## Lettura dei dati sperimentali da Flubacher 1960

        Non è corretto confrontare direttamente i dati di Flubacher, perché sono capacità termiche a pressione costante. Io cerco invece quelle a volume costante. Queste sono disponibili nell'articolo di Holzapfel e Klotz del 2021.
        """
    )
    return


@app.cell
def __():
    import numpy as np

    flubacher = np.loadtxt("../flubacher.txt")
    return flubacher, np


@app.cell
def __(flubacher):
    flubacher_temps = flubacher[:, 0]
    flubacher_heat_capacities_cal = flubacher[:, 1]
    cal2J = 4.184
    flubacher_heat_capacities_J = flubacher_heat_capacities_cal * cal2J
    return (
        cal2J,
        flubacher_heat_capacities_J,
        flubacher_heat_capacities_cal,
        flubacher_temps,
    )


@app.cell
def __(
    flubacher_heat_capacities_J,
    flubacher_temps,
    plt,
    thermal_properties,
):
    plt.scatter(
        x=flubacher_temps, y=flubacher_heat_capacities_J, label="Flubacher 1960"
    )
    plt.scatter(
        x=thermal_properties["temperatures"],
        y=thermal_properties["heat_capacity"],
        label="MACE-MP-0",
    )
    plt.xlabel("Temperature (K)")
    plt.ylabel("Heat capacity (J/K/mol)")

    plt.xlim(0, 50)

    plt.legend()
    return


@app.cell
def __(df, flubacher_heat_capacities_J, flubacher_temps, plt):
    x = flubacher_temps**2
    y = flubacher_heat_capacities_J / flubacher_temps**3 * 1e5

    plt.scatter(x, y)

    plt.scatter(
        df["temperatures"] ** 2,
        df["heat_capacity"] / df["temperatures"] ** 3 * 1e5,
    )

    plt.xlabel("$T^2$")
    plt.ylabel("$10^5 C/T^3$")

    plt.xlim(-100,2000)
    plt.gca()
    return x, y


@app.cell
def __(thermal_properties):
    thermal_properties
    return


@app.cell
def __(thermal_properties):
    import pandas as pd
    df = pd.DataFrame(thermal_properties).dropna()
    df
    return df, pd


@app.cell
def __(holzapfel, mo):
    mo.hstack(
        [
            mo.md(
                """
        ## Lettura dei dati di Holzapfel e Klotz 2021

        Presi dalla TABLE V.

        
        """
            ),
            holzapfel,
        ]
    )
    return


@app.cell
def __(pd):
    holzapfel = pd.read_csv(
        "../holzapfel_klotz_2021.txt", delimiter="[\ \t]+", engine="python"
    )
    return holzapfel,


@app.cell
def __(R_KJmol, holzapfel, plt, thermal_properties):
    plt.scatter(x=holzapfel["T(K)"], y=holzapfel["C_V/R"], label="Holzapfel 2021")
    plt.scatter(
        x=thermal_properties["temperatures"],
        y=thermal_properties["heat_capacity"] / R_KJmol,
        label="MACE-MP-0",
    )
    plt.legend()
    return


@app.cell
def __(R_KJmol, holzapfel, np, plt, thermal_properties):
    plt.scatter(
        x=np.log10(holzapfel["T(K)"].loc[1:]),
        y=np.log10(holzapfel["C_V/R"].loc[1:]),
        label="Holzapfel",
    )
    plt.scatter(
        x=np.log10(thermal_properties["temperatures"][1:]),
        y=np.log10(thermal_properties["heat_capacity"][1:] / R_KJmol),
        label="MACE-MP-0",
    )
    plt.legend()
    return


@app.cell
def __(np):
    help(np.log)
    return


if __name__ == "__main__":
    app.run()
