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
        start=0, stop=999, debounce=False, label="t_min", value=0, show_value=True
    ).form()
    t_max = mo.ui.slider(
        start=1,
        stop=1000,
        debounce=False,
        label="t_max",
        value=100,
        show_value=True,
    ).form()
    t_step = mo.ui.slider(
        start=1,
        stop=100,
        debounce=False,
        label="t_step",
        value=10,
        show_value=True,
    ).form()
    return t_max, t_min, t_step


@app.cell
def __(mo, t_max, t_min, t_step):
    mo.hstack([t_min, t_max, t_step])
    return


@app.cell
def __(ph, plt, t_max, t_min, t_step):
    ph.run_mesh()
    ph.run_thermal_properties(
        t_step=t_step.value, t_min=t_min.value, t_max=t_max.value
    )
    thermal_properties = ph.get_thermal_properties_dict()

    plt.plot(
        thermal_properties["temperatures"],
        thermal_properties["heat_capacity"],
        marker="o",
    )

    plt.xlabel("Temperature (K)")
    # Unità di misura definite in
    # https://phonopy.github.io/phonopy/setting-tags.html#thermal-properties-related-tags
    plt.ylabel("Heat capacity (J/K/mol)")
    return thermal_properties,


if __name__ == "__main__":
    app.run()
