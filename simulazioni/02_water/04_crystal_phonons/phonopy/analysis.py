import marimo

__generated_with = "0.6.22"
app = marimo.App()


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        """
        # Studio dei modi vibrazionali del polimorfo Ih del ghiaccio

        Si sono calcolate le forze e la matrice delle costanti di forza col metodo degli spostamenti finiti, usando una supercella \(2 \\times 2 \\times 2\) e \(3\\times3\\times3\).
        Queste quantità fisiche sono conservate nel file `phonopy_params.yaml`.

        Con gli strumenti di ASE è possibile visualizzare la prima zona di Brillouin della geometria in esame, ottenendo le coordinate dei punti speciali.

        Si confrontano i risultati delle simulazioni con i valori in letteratura disponibili.
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

    from mpl_toolkits.axes_grid1 import ImageGrid
    from phonopy.phonon.band_structure import BandPlot
    from phonopy.phonon.band_structure import BandStructure
    from phonopy.phonon.band_structure import band_plot
    return (
        BandPlot,
        BandStructure,
        ImageGrid,
        band_plot,
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
def __(atoms, mo):
    mo.md(
        f"""
        **atoms.cell.reciprocal**:  
        {atoms.cell.reciprocal()}
        """
    )
    return


@app.cell
def __(bandpath, mo):
    mo.md(
        f"""
        **Bandpath**:  
        {bandpath}

        **Special points**:
        """
    )
    return


@app.cell
def __(bandpath):
    bandpath.plot()
    return


@app.cell
def __(atoms, bandpath, mo):
    _labels = []
    _path = [bandpath.special_points.get(key) for key in _labels]

    _bandpath = atoms.cell.bandpath(_labels)

    mo.mpl.interactive(_bandpath.plot())
    return


@app.cell
def __(mo):
    mo.md(
        """
        ## Costruzione di Strässle

        Bisogna costruire un percorso. Scegliamo per ora di seguire quello studiato nell'articolo di Strässle, Saitta, Klotz del 2004, che segue il percorso \(A\Gamma K M \Gamma\).
        """
    )
    return


@app.cell
def __(atoms, bandpath, get_band_qpoints_and_path_connections, mo):
    labels = ["A", "G", "K", "M", "G"]
    path = [bandpath.special_points.get(key) for key in labels]

    _bandpath = atoms.cell.bandpath(labels)

    qpoints, connections = get_band_qpoints_and_path_connections(
        # qpoints = get_band_qpoints(
        [path],
        npoints=101,
        # rec_lattice=atoms.cell.reciprocal()
    )

    mo.mpl.interactive(_bandpath.plot())
    return connections, labels, path, qpoints


@app.cell
def __(mo):
    mo.md(
        rf"""
        ## Costruzione di Gupta

        L'articolo @guptaPhononsAnomalousThermal2018
        """
    )
    return


@app.cell
def __(atoms, mo):
    labels_gupta = ["G", "A", "K", "H", "M", "L", "G"]
    bandpath_gupta = atoms.cell.bandpath(labels_gupta)
    mo.mpl.interactive(bandpath_gupta.plot())
    return bandpath_gupta, labels_gupta


@app.cell
def __(
    bandpath_gupta,
    get_band_qpoints_and_path_connections,
    labels_gupta,
):
    # Build the path for phonopy
    path_gupta = [bandpath_gupta.special_points.get(key) for key in labels_gupta]

    qpoints_gupta, connections_gupta = get_band_qpoints_and_path_connections(
        [path_gupta],
        npoints=101,
    )

    path_gupta
    return connections_gupta, path_gupta, qpoints_gupta


@app.cell
def __(THz2meV, mace_ice13_1_s3, meV2THz, mo, os, plt):
    _fig, _ax = plt.subplots(
        ncols=2, nrows=1, layout="constrained", width_ratios=[9999, 1]
    )
    mace_ice13_1_s3["ph"].band_structure.plot(_ax)

    _ax[0].set_ylim(bottom=0)
    _ax[0].grid(axis="x")
    _ax[0].set_ylabel("Frequency (THz)")

    _fig.delaxes(_ax[1])
    plt.title("MACE-ICE13-1 supercell 3x3x3")

    # Remove top frame
    _ax[0].spines["top"].set_visible(False)

    _secax = _ax[0].secondary_yaxis("right", functions=(THz2meV, meV2THz))
    _secax.set_ylabel("E (meV)")

    os.makedirs("Grafici", exist_ok=True)
    # plt.savefig("Grafici/bandstructure_mace-ice13-1_s3_gupta_full.svg")

    mo.mpl.interactive(_fig)
    return


@app.cell
def __(THz2meV, mace_ice13_1_s3, meV2THz, mo, plt):
    _fig, _ax = plt.subplots(
        ncols=2, nrows=1, layout="constrained", width_ratios=[9999, 1]
    )
    mace_ice13_1_s3["ph"].band_structure.plot(_ax)

    _ax[0].set_ylim(0, 5)
    _ax[0].grid(axis="x")
    _ax[0].set_ylabel("Frequency (THz)")


    _fig.delaxes(_ax[1])
    plt.title("MACE-ICE13-1 supercell 3x3x3")

    # Remove top frame
    _ax[0].spines["top"].set_visible(False)
    # _ax[0].spines["left"].set_visible(False)
    # _ax[0].spines["right"].set_visible(False)

    _secax = _ax[0].secondary_yaxis("right", functions=(THz2meV, meV2THz))
    _secax.set_ylabel("E (meV)")

    import os

    os.makedirs("Grafici", exist_ok=True)
    # plt.savefig("Grafici/bandstructure_mace-ice13-1_s3_gupta.svg")

    mo.mpl.interactive(_fig)
    return os,


@app.cell
def __(THz2meV, band_plot, meV2THz):
    def plot_the_structure(ax, ph, connections, labels, fmt):
        _dict = ph.get_band_structure_dict()
        _frequencies = _dict["frequencies"]
        _distances = _dict["distances"]
        band_plot(ax, _frequencies, _distances, connections, labels, fmt=fmt)

    def style_bandstructure_plot(ax):
        ax.set_ylim(0, 5)
        ax.grid(axis="x")
        ax.set_ylabel("Frequency (THz)")
        ax.spines["top"].set_visible(False)
        
    def make_second_axis(ax):
        _secax = ax.secondary_yaxis("right", functions=(THz2meV, meV2THz))
        _secax.set_ylabel("E (meV)")
    return make_second_axis, plot_the_structure, style_bandstructure_plot


@app.cell
def __(
    Line2D,
    connections_gupta,
    labels_gupta,
    mace_ice13_1,
    mace_ice13_1_s3,
    make_second_axis,
    mo,
    plot_the_structure,
    plt,
    style_bandstructure_plot,
):
    # mace_ice13_1["ph"].run_band_structure(
    #     qpoints_gupta,
    #     path_connections=connections_gupta,
    #     labels=labels_gupta,
    #     is_legacy_plot=False,
    # )

    _fig, _axs = plt.subplots(
        ncols=2, nrows=1, layout="constrained", width_ratios=[9999, 1]
    )


    plot_the_structure(
        _axs, mace_ice13_1["ph"], connections_gupta, labels_gupta, fmt="r--"
    )
    plot_the_structure(
        _axs, mace_ice13_1_s3["ph"], connections_gupta, labels_gupta, fmt="g-"
    )

    _fig.delaxes(_axs[1])


    style_bandstructure_plot(_axs[0])

    _custom_lines = [
        Line2D([0], [0], color="r", lw=2, ls="--"),
        Line2D([0], [0], color="g", lw=2),
    ]
    _axs[0].legend(_custom_lines, ["MACE-ICE13-1 S2", "MACE-ICE13-1 S3"])


    make_second_axis(_axs[0])
    # _secax = _axs[0].secondary_yaxis("right", functions=(THz2meV, meV2THz))
    # _secax.set_ylabel("E (meV)")

    # _fig.savefig("Grafici/bandstructure_mace-ice13-1_s2-s3_gupta_zoom.svg")
    mo.mpl.interactive(_fig)
    return


@app.cell
def __(
    Line2D,
    connections_gupta,
    labels_gupta,
    mace_ice13_1,
    mace_mp_0,
    make_second_axis,
    mo,
    plot_the_structure,
    plt,
    style_bandstructure_plot,
):
    # mace_mp_0["ph"].run_band_structure(
    #     qpoints_gupta,
    #     path_connections=connections_gupta,
    #     labels=labels_gupta,
    #     is_legacy_plot=False,
    # )

    _fig, _axs = plt.subplots(
        ncols=2, nrows=1, layout="constrained", width_ratios=[9999, 1]
    )

    plot_the_structure(
        _axs, mace_ice13_1["ph"], connections_gupta, labels_gupta, fmt="r-"
    )
    plot_the_structure(
        _axs, mace_mp_0["ph"], connections_gupta, labels_gupta, fmt="b-"
    )

    _fig.delaxes(_axs[1])

    style_bandstructure_plot(_axs[0])

    _custom_lines = [
        Line2D([0], [0], color="r", lw=2),
        Line2D([0], [0], color="b", lw=2),
    ]
    _axs[0].legend(_custom_lines, ["MACE-ICE13-1", "MACE-MP-0 medium + D3"])

    make_second_axis(_axs[0])

    _fig.suptitle("Supercell $2 \\times 2 \\times 2$")

    # _fig.savefig("Grafici/bandstructure_mace-ice13-1-mace-mp.0_s2_gupta_zoom.svg")
    mo.mpl.interactive(_fig)
    return


@app.cell
def __(mace_mp_0, make_second_axis, mo, os, plt, style_bandstructure_plot):
    _fig, _ax = plt.subplots(
        ncols=2, nrows=1, layout="constrained", width_ratios=[9999, 1]
    )
    mace_mp_0["ph"].band_structure.plot(_ax)

    _fig.delaxes(_ax[1])

    style_bandstructure_plot(_ax[0])

    plt.title("MACE-MP-0 supercell 2x2x2")

    make_second_axis(_ax[0])

    os.makedirs("Grafici", exist_ok=True)
    # plt.savefig("Grafici/bandstructure_mace-mp-0_s2_gupta.svg")

    mo.mpl.interactive(_fig)
    return


@app.cell
def __():
    return


@app.cell
def __(connections, labels, phonopy, qpoints):
    mace_mp_0 = {}
    mace_mp_0["basepath"] = "MACE-MP-0-d3-medium/ICE-Ih/supercell=2/"
    mace_mp_0["ph"] = phonopy.load(
        mace_mp_0["basepath"] + "phonopy_params.yaml"
    )
    mace_mp_0["ph"].run_band_structure(
        qpoints, path_connections=connections, labels=labels, is_legacy_plot=False
    )
    return mace_mp_0,


@app.cell
def __(connections, labels, phonopy, qpoints):
    mace_ice13_1 = {}
    mace_ice13_1["basepath"] = "MACE-ICE13-1/ICE-Ih/supercell=2/"
    mace_ice13_1["ph"] = phonopy.load(
        mace_ice13_1["basepath"] + "phonopy_params.yaml"
    )
    mace_ice13_1["ph"].run_band_structure(
        qpoints, path_connections=connections, labels=labels, is_legacy_plot=False
    )
    return mace_ice13_1,


@app.cell
def __():
    # enhance below with https://docs.marimo.io/guides/performance.html
    import functools
    return functools,


@app.cell
def __(functools, phonopy):
    # WARNING Expensive, run just once
    mace_ice13_1_s3 = {}
    mace_ice13_1_s3["basepath"] = "MACE-ICE13-1/ICE-Ih/supercell=3/"

    @functools.cache
    def phonopy_load_cached(input):
        return phonopy.load(input)

    mace_ice13_1_s3["ph"] = phonopy_load_cached(
        mace_ice13_1_s3["basepath"] + "phonopy_params.yaml"
    )
    return mace_ice13_1_s3, phonopy_load_cached


@app.cell
def __(mo):
    mo.md("## Plot di tutte le bande trovate")
    return


@app.cell
def __(connections, labels, mace_ice13_1_s3, qpoints):
    mace_ice13_1_s3["ph"].run_band_structure(
        qpoints, path_connections=connections, labels=labels, is_legacy_plot=False
    )
    mace_ice13_1_s3["ph"].plot_band_structure().gca()
    # plt.savefig("mace_ice13_1_s3_band_structure_complete.svg")
    return


@app.cell
def __(connections_gupta, labels_gupta, mace_ice13_1_s3, qpoints_gupta):
    mace_ice13_1_s3["ph"].run_band_structure(
        qpoints_gupta,
        path_connections=connections_gupta,
        labels=labels_gupta,
        is_legacy_plot=False,
    )
    mace_ice13_1_s3["ph"].plot_band_structure().gca()
    return


@app.cell
def __(mo):
    mo.md("## Zoom sulle bande più basse")
    return


@app.cell
def __(mace_ice13_1, plt):
    _fig, _ax = plt.subplots(
        ncols=2, nrows=1, layout="constrained", width_ratios=[9999, 1]
    )
    mace_ice13_1["ph"].band_structure.plot(_ax)

    _ax[0].set_ylim(0, 5)
    _ax[0].grid(axis="x")
    _ax[0].set_ylabel("Frequency (THz)")

    _fig.delaxes(_ax[1])
    _fig
    return


@app.cell
def __(mace_ice13_1_s3, plt):
    _fig, _ax = plt.subplots(
        ncols=2, nrows=1, layout="constrained", width_ratios=[9999, 1]
    )
    mace_ice13_1_s3["ph"].band_structure.plot(_ax)

    _ax[0].set_ylim(0, 5)
    _ax[0].grid(axis="x")
    _ax[0].set_ylabel("Frequency (THz)")

    _fig.delaxes(_ax[1])
    plt.title("MACE-ICE13-1 supercell 3x3x3")
    _fig
    return


@app.cell
def __(band_plot, connections, labels, mace_ice13_1, mace_mp_0, plt):
    n = len(
        [x for x in mace_ice13_1["ph"]._band_structure.path_connections if not x]
    )

    _fig, _axs = plt.subplots(
        ncols=2, nrows=1, layout="constrained", width_ratios=[9999, 1]
    )
    # _fig = plt.figure(layout="constrained")
    # _axs = ImageGrid(_fig, 111, nrows_ncols=(1, n), axes_pad=0.11, label_mode="L")
    # mace_ice13_1["ph"]._band_structure.plot(_axs)
    _dict = mace_ice13_1["ph"].get_band_structure_dict()
    _frequencies = _dict["frequencies"]
    _distances = _dict["distances"]

    band_plot(_axs, _frequencies, _distances, connections, labels, fmt="r-")

    _dict = mace_mp_0["ph"].get_band_structure_dict()
    _frequencies = _dict["frequencies"]
    _distances = _dict["distances"]
    band_plot(_axs, _frequencies, _distances, connections, labels, fmt="b-")

    _axs[0].set_ylim(0, 5)
    _fig.delaxes(_axs[1])
    _axs[0].set_ylabel("Frequency (THz)")
    _axs[0].grid()

    from matplotlib.lines import Line2D
    custom_lines = [
        Line2D([0], [0], color="r", lw=2),
        Line2D([0], [0], color="b", lw=2),
    ]
    _axs[0].legend(custom_lines, ["MACE-ICE13-1", "MACE-MP-0"])

    _fig.suptitle("Supercell $2 \\times 2 \\times 2$")

    _fig
    return Line2D, custom_lines, n


@app.cell
def __(
    Line2D,
    band_plot,
    connections,
    labels,
    mace_ice13_1,
    mace_ice13_1_s3,
    plt,
):
    _n = len(
        [x for x in mace_ice13_1["ph"]._band_structure.path_connections if not x]
    )

    _fig, _axs = plt.subplots(
        ncols=2, nrows=1, layout="constrained", width_ratios=[9999, 1]
    )
    _dict = mace_ice13_1["ph"].get_band_structure_dict()
    _frequencies = _dict["frequencies"]
    _distances = _dict["distances"]

    band_plot(_axs, _frequencies, _distances, connections, labels, fmt="r--")

    _dict = mace_ice13_1_s3["ph"].get_band_structure_dict()
    _frequencies = _dict["frequencies"]
    _distances = _dict["distances"]
    band_plot(_axs, _frequencies, _distances, connections, labels, fmt="g-")

    _axs[0].set_ylim(0, 5)
    _fig.delaxes(_axs[1])
    _axs[0].set_ylabel("Frequency (THz)")
    _axs[0].grid()

    _custom_lines = [
        Line2D([0], [0], color="r", lw=2, ls="--"),
        Line2D([0], [0], color="g", lw=2),
    ]
    _axs[0].legend(_custom_lines, ["MACE-ICE13-1 S2", "MACE-ICE13-1 S3"])

    # _fig.savefig("mace_ice13_1_s2vss3_band_structure_zoom.svg")
    _fig
    return


@app.cell
def __(mo):
    mesh = mo.ui.slider(1, 16, value=4, show_value=True).form()
    return mesh,


@app.cell
def __(mesh, mo):
    mo.md(
        f"""
        ## Calcolo della DOS

        Scegli di seguito il lato della griglia per il calcolo della DOS:

        {mesh}
        """
    )
    return


@app.cell
def __(mace_ice13_1, mace_ice13_1_s3, mace_mp_0, mesh):
    _sigma = 0.05

    mace_ice13_1["ph"].run_mesh([mesh.value, mesh.value, mesh.value])
    mace_ice13_1["ph"].run_total_dos(
        sigma=_sigma,
        freq_min=0,
        freq_max=14,
        freq_pitch=None,
        use_tetrahedron_method=True,
    )
    # ph.write_total_dos(filename=f"total_dos_mesh_{mesh.value}.dat")

    mace_ice13_1_s3["ph"].run_mesh([mesh.value, mesh.value, mesh.value])
    mace_ice13_1_s3["ph"].run_total_dos(
        sigma=_sigma,
        freq_min=0,
        freq_max=14,
        freq_pitch=None,
        use_tetrahedron_method=True,
    )

    mace_mp_0["ph"].run_mesh([mesh.value, mesh.value, mesh.value])
    mace_mp_0["ph"].run_total_dos(
        sigma=_sigma,
        freq_min=0,
        freq_max=14,
        freq_pitch=None,
        use_tetrahedron_method=True,
    )
    return


@app.cell
def __(mace_ice13_1_s3, mace_mp_0, mesh, plt):
    fig, ax = plt.subplots(layout="constrained")

    # ax.plot(
    #     mace_ice13_1["ph"].get_total_dos_dict()["frequency_points"],
    #     mace_ice13_1["ph"].get_total_dos_dict()["total_dos"],
    #     color="red",
    #     label="MACE-ICE13-1",
    #     zorder=10,
    # )

    ax.plot(
        mace_ice13_1_s3["ph"].get_total_dos_dict()["frequency_points"],
        mace_ice13_1_s3["ph"].get_total_dos_dict()["total_dos"],
        color="green",
        label="MACE-ICE13-1 S=3",
        zorder=20,
    )

    ax.plot(
        mace_mp_0["ph"].get_total_dos_dict()["frequency_points"],
        mace_mp_0["ph"].get_total_dos_dict()["total_dos"],
        color="blue",
        label="MACE-MP-0",
        zorder=1,
        linestyle="-",
        linewidth=0.5,
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
    plt.grid()
    plt.legend(loc="upper left")
    ax
    return K2THz, THz2K, ax, fig, secax


@app.cell
def __(mo):
    mo.md(
        """
        ## Calcolo dell'energia in meV per confronto con _delrosso2021_

        La Figure 1 dell'articolo arriva orientativamente a 130meV.
        Questo corrisponde a:

        ```
        ❯ pint-convert 130meV THz
        130 millielectron_volt = 31.4338601471 THz
        ❯ pint-convert meV THz   
        1 millielectron_volt = 0.241798924208 THz
        ❯ pint-convert THz meV
        1 terahertz = 4.13566769692 meV
        ```

        Bisogna quindi calcolare la DOS fino a quella frequenza massima.
        """
    )
    return


@app.cell
def __(mace_ice13_1_s3, plt):
    _mesh = 32
    mace_ice13_1_s3["ph"].run_mesh([_mesh] * 3)

    mace_ice13_1_s3["ph"].run_total_dos(
        sigma=0.05,
        freq_min=0,
        freq_max=36,
        freq_pitch=None,
        use_tetrahedron_method=False,
    )

    _fig, _ax = plt.subplots(layout="constrained")

    _ax.plot(
        mace_ice13_1_s3["ph"].get_total_dos_dict()["frequency_points"],
        mace_ice13_1_s3["ph"].get_total_dos_dict()["total_dos"],
        color="green",
        label="MACE-ICE13-1 S=3",
        zorder=20,
    )

    _ax.set_ylim(bottom=0)
    _ax.set_xlim(left=0)
    _ax.set_xlabel("Freq (THz)")
    _ax.set_ylabel("DOS (1/THz)")
    _ax.set_title(f"{_mesh}x{_mesh}x{_mesh} mesh")


    def THz2meV(THz):
        return THz * 4.13566769692


    def meV2THz(meV):
        return meV * 0.241798924208


    _secax = _ax.secondary_xaxis("top", functions=(THz2meV, meV2THz))
    _secax.set_xlabel("E (meV)")
    plt.grid()
    plt.legend()

    # _fig.savefig("mace_ice13_1_s3_dos.svg")

    _ax
    return THz2meV, meV2THz


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
def __(mace_ice13_1, mace_ice13_1_s3, mace_mp_0, numero_molecole, plt):
    mace_ice13_1["ph"].run_mesh(mesh=100.0)
    # ph.run_mesh(mesh=[_m] * 3)


    def build_thermal_properties(ph, t_step=1, t_min=0, t_max=50):
        ph.run_thermal_properties(t_step=t_step, t_min=t_min, t_max=t_max)
        tp = ph.get_thermal_properties_dict()
        tp["heat_capacity"] = tp["heat_capacity"] / numero_molecole
        return tp


    mace_ice13_1["thermal_properties_flubacher"] = build_thermal_properties(
        mace_ice13_1["ph"]
    )
    mace_mp_0["thermal_properties_flubacher"] = build_thermal_properties(
        mace_mp_0["ph"]
    )
    mace_ice13_1_s3["thermal_properties_flubacher"] = build_thermal_properties(
        mace_ice13_1_s3["ph"]
    )

    R_KJmol = 8.31446261815

    plt.plot(
        mace_ice13_1["thermal_properties_flubacher"]["temperatures"],
        mace_ice13_1["thermal_properties_flubacher"]["heat_capacity"] / R_KJmol,
        marker="o",
        label="MACE-ICE13-1",
        color="red",
    )

    plt.plot(
        mace_mp_0["thermal_properties_flubacher"]["temperatures"],
        mace_mp_0["thermal_properties_flubacher"]["heat_capacity"] / R_KJmol,
        marker="o",
        label="MACE-MP-0",
        color="blue",
    )

    plt.plot(
        mace_ice13_1_s3["thermal_properties_flubacher"]["temperatures"],
        mace_ice13_1_s3["thermal_properties_flubacher"]["heat_capacity"] / R_KJmol,
        marker="o",
        label="MACE-ICE13-1 S3",
        color="green",
    )

    plt.xlabel("Temperature (K)")
    # Unità di misura definite in
    # https://phonopy.github.io/phonopy/setting-tags.html#thermal-properties-related-tags
    # you have to divide the value by number of formula unit in your unit cell by yourself
    plt.ylabel("Heat capacity (J/K/mol)")
    plt.legend()
    return R_KJmol, build_thermal_properties


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

    flubacher = np.loadtxt("analysis_assets/flubacher.txt")
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
    mace_ice13_1,
    mace_ice13_1_s3,
    mace_mp_0,
    plt,
):
    plt.scatter(
        x=flubacher_temps,
        y=flubacher_heat_capacities_J,
        label="Flubacher 1960",
        color="green",
        marker="x"
    )
    plt.scatter(
        x=mace_ice13_1["thermal_properties_flubacher"]["temperatures"],
        y=mace_ice13_1["thermal_properties_flubacher"]["heat_capacity"],
        label="MACE-ICE13-1",
        color="red",
    )

    plt.scatter(
        x=mace_mp_0["thermal_properties_flubacher"]["temperatures"],
        y=mace_mp_0["thermal_properties_flubacher"]["heat_capacity"],
        label="MACE-MP-0",
        color="blue",
    )

    plt.scatter(
        x=mace_ice13_1_s3["thermal_properties_flubacher"]["temperatures"],
        y=mace_ice13_1_s3["thermal_properties_flubacher"]["heat_capacity"],
        label="MACE-ICE13-1 S3",
        color="green",
    )

    plt.xlabel("Temperature (K)")
    plt.ylabel("Heat capacity (J/K/mol)")

    plt.xlim(-0.5, 31)
    plt.ylim(-0.2, 4)

    plt.legend()
    return


@app.cell(hide_code=True)
def __(
    R_KJmol,
    flubacher_heat_capacities_J,
    flubacher_temps,
    holzapfel,
    mace_ice13_1,
    mace_ice13_1_s3,
    mace_mp_0,
    plt,
):
    x = flubacher_temps**2
    y = flubacher_heat_capacities_J / flubacher_temps**3 * 1e5

    plt.scatter(x, y, color="green", label="Flubacher 1960", marker="x")

    plt.scatter(
        x=holzapfel["T(K)"] ** 2,
        y=holzapfel["C_V/R"] * R_KJmol / holzapfel["T(K)"] ** 3 * 1e5,
        label="Holzapfel 2021",
        color="green",
        marker="+",
    )

    plt.scatter(
        mace_ice13_1_s3["thermal_properties_flubacher"]["temperatures"] ** 2,
        mace_ice13_1_s3["thermal_properties_flubacher"]["heat_capacity"]
        / mace_ice13_1_s3["thermal_properties_flubacher"]["temperatures"] ** 3
        * 1e5,
        color="green",
        label="MACE-ICE13-1 S=3",
        marker="s",
    )

    plt.scatter(
        mace_ice13_1["thermal_properties_flubacher"]["temperatures"] ** 2,
        mace_ice13_1["thermal_properties_flubacher"]["heat_capacity"]
        / mace_ice13_1["thermal_properties_flubacher"]["temperatures"] ** 3
        * 1e5,
        color="red",
        label="MACE-ICE13-1 S=2",
    )

    plt.scatter(
        mace_mp_0["thermal_properties_flubacher"]["temperatures"] ** 2,
        mace_mp_0["thermal_properties_flubacher"]["heat_capacity"]
        / mace_mp_0["thermal_properties_flubacher"]["temperatures"] ** 3
        * 1e5,
        color="blue",
        label="MACE-MP-0",
    )

    plt.xlabel("$T^2 \, \mathrm{(K^2)}$")
    plt.ylabel("$10^5 C/T^3 \, \mathrm{(J/K^4/mol)}$")

    plt.xlim(-10, 800)
    plt.legend()
    plt.grid()
    plt.title("Ice Ih heat capacity")

    # plt.savefig("heat_capacity_t3.svg")
    plt.gca()
    return x, y


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
def __():
    import pandas as pd

    holzapfel = pd.read_csv(
        "analysis_assets/holzapfel_klotz_2021.txt",
        delimiter="[\ \t]+",
        engine="python",
    )
    return holzapfel, pd


@app.cell
def __(
    R_KJmol,
    build_thermal_properties,
    holzapfel,
    mace_ice13_1,
    mace_ice13_1_s3,
    mace_mp_0,
    plt,
):
    plt.scatter(
        x=holzapfel["T(K)"],
        y=holzapfel["C_V/R"],
        label="Holzapfel 2021",
        color="green",
        marker="x",
    )

    mace_ice13_1["thermal_properties"] = build_thermal_properties(
        mace_ice13_1["ph"], t_step=5, t_min=0, t_max=265
    )


    mace_ice13_1_s3["thermal_properties"] = build_thermal_properties(
        mace_ice13_1_s3["ph"], t_step=5, t_min=0, t_max=265
    )

    mace_mp_0["thermal_properties"] = build_thermal_properties(
        mace_mp_0["ph"], t_step=5, t_min=0, t_max=265
    )

    plt.scatter(
        x=mace_ice13_1_s3["thermal_properties"]["temperatures"],
        y=mace_ice13_1_s3["thermal_properties"]["heat_capacity"] / R_KJmol,
        label="MACE-ICE13-1 S3",
        color="green",
        marker="s",
        s=25,
    )

    plt.scatter(
        x=mace_ice13_1["thermal_properties"]["temperatures"],
        y=mace_ice13_1["thermal_properties"]["heat_capacity"] / R_KJmol,
        label="MACE-ICE13-1",
        color="red",
        marker="o",
        s=25,
    )

    plt.scatter(
        x=mace_mp_0["thermal_properties"]["temperatures"],
        y=mace_mp_0["thermal_properties"]["heat_capacity"] / R_KJmol,
        label="MACE-MP-0",
        color="blue",
        marker="o",
        s=25,
    )

    plt.xlabel("Temperature (K)")
    plt.ylabel("$C_V/R$")
    plt.title("Ice Ih heat capacity")
    plt.legend()
    plt.grid()
    # plt.savefig("heat_capacity_all_temps.svg")
    plt.gca()
    return


@app.cell(hide_code=True)
def __(
    R_KJmol,
    build_thermal_properties,
    holzapfel,
    mace_ice13_1,
    mace_ice13_1_s3,
    mace_mp_0,
    plt,
):
    # Zoom sulla zona a basse temperature

    plt.scatter(
        x=holzapfel["T(K)"],
        y=holzapfel["C_V/R"],
        label="Holzapfel 2021",
        color="green",
        marker="x",
    )

    mace_ice13_1["thermal_properties"] = build_thermal_properties(
        mace_ice13_1["ph"], t_step=5, t_min=0, t_max=265
    )

    plt.scatter(
        x=mace_ice13_1["thermal_properties"]["temperatures"],
        y=mace_ice13_1["thermal_properties"]["heat_capacity"] / R_KJmol,
        label="MACE-ICE13-1",
        color="red",
    )

    mace_mp_0["thermal_properties"] = build_thermal_properties(
        mace_mp_0["ph"], t_step=5, t_min=0, t_max=265
    )

    plt.scatter(
        x=mace_mp_0["thermal_properties"]["temperatures"],
        y=mace_mp_0["thermal_properties"]["heat_capacity"] / R_KJmol,
        label="MACE-MP-0",
        color="blue",
    )

    mace_ice13_1_s3["thermal_properties"] = build_thermal_properties(
        mace_ice13_1_s3["ph"], t_step=5, t_min=0, t_max=265
    )

    plt.scatter(
        x=mace_ice13_1_s3["thermal_properties"]["temperatures"],
        y=mace_ice13_1_s3["thermal_properties"]["heat_capacity"] / R_KJmol,
        label="MACE-ICE13-1 S3",
        color="green",
    )

    plt.xlabel("Temperature (K)")
    plt.ylabel("$C_V/R$")
    plt.title("Constant volume heat capacity")

    plt.xlim(-0.5, 31)
    plt.ylim(-0.01, 0.5)

    plt.legend()
    return


@app.cell
def __(
    R_KJmol,
    build_thermal_properties,
    flubacher_heat_capacities_J,
    flubacher_temps,
    holzapfel,
    mace_ice13_1,
    mace_ice13_1_s3,
    mace_mp_0,
    plt,
):
    # Plot di Holzapfel insieme a Flubacher

    mace_ice13_1_s3["thermal_properties"] = build_thermal_properties(
        mace_ice13_1_s3["ph"], t_step=1, t_min=0, t_max=30
    )
    mace_ice13_1["thermal_properties"] = build_thermal_properties(
        mace_ice13_1["ph"], t_step=1, t_min=0, t_max=30
    )
    mace_mp_0["thermal_properties"] = build_thermal_properties(
        mace_mp_0["ph"], t_step=1, t_min=0, t_max=30
    )

    plt.scatter(
        x=holzapfel["T(K)"],
        y=holzapfel["C_V/R"] * R_KJmol,
        label="Holzapfel 2021",
        color="green",
        marker="+",
    )

    plt.scatter(
        x=flubacher_temps,
        y=flubacher_heat_capacities_J,
        label="Flubacher 1960",
        color="green",
        marker="x",
        s=25,
    )

    plt.plot(
        mace_ice13_1_s3["thermal_properties"]["temperatures"],
        mace_ice13_1_s3["thermal_properties"]["heat_capacity"],
        label="MACE-ICE13-1 S=3",
        color="green",
        marker="s",
        markersize=4,
    )

    plt.plot(
        mace_ice13_1["thermal_properties"]["temperatures"],
        mace_ice13_1["thermal_properties"]["heat_capacity"],
        label="MACE-ICE13-1 S=2",
        color="red",
        linestyle="--",
        marker="o",
        markersize=4,
    )

    plt.plot(
        mace_mp_0["thermal_properties"]["temperatures"],
        mace_mp_0["thermal_properties"]["heat_capacity"],
        label="MACE-MP-0",
        color="blue",
        marker="o",
        markersize=4,
    )

    plt.xlabel("Temperature (K)")
    plt.ylabel("J/K/mol")
    plt.title("Ice Ih heat capacity")

    plt.xlim(0, 31)
    plt.ylim(0, 4)

    plt.legend()
    # plt.savefig("heat_capacity_small_temp.svg")
    plt.gca()
    return


@app.cell
def __(
    R_KJmol,
    build_thermal_properties,
    holzapfel,
    mace_ice13_1,
    mace_ice13_1_s3,
    mace_mp_0,
    np,
    plt,
):
    mace_ice13_1_s3["thermal_properties"] = build_thermal_properties(
        mace_ice13_1_s3["ph"], t_step=1, t_min=0, t_max=250
    )
    mace_ice13_1["thermal_properties"] = build_thermal_properties(
        mace_ice13_1["ph"], t_step=1, t_min=0, t_max=250
    )
    mace_mp_0["thermal_properties"] = build_thermal_properties(
        mace_mp_0["ph"], t_step=1, t_min=0, t_max=250
    )

    plt.scatter(
        x=np.log10(holzapfel["T(K)"].loc[1:]),
        y=np.log10(holzapfel["C_V/R"].loc[1:]),
        label="Holzapfel 2021",
        color="green",
        marker="x",
    )

    plt.scatter(
        x=np.log10(mace_ice13_1_s3["thermal_properties"]["temperatures"][1:]),
        y=np.log10(
            mace_ice13_1_s3["thermal_properties"]["heat_capacity"][1:] / R_KJmol
        ),
        label="MACE-ICE13-1 S3",
        color="green",
        s=10,
        marker="s",
    )

    plt.scatter(
        x=np.log10(mace_ice13_1["thermal_properties"]["temperatures"][1:]),
        y=np.log10(
            mace_ice13_1["thermal_properties"]["heat_capacity"][1:] / R_KJmol
        ),
        label="MACE-ICE13-1",
        color="red",
        s=10,
    )

    plt.scatter(
        x=np.log10(mace_mp_0["thermal_properties"]["temperatures"][1:]),
        y=np.log10(mace_mp_0["thermal_properties"]["heat_capacity"][1:] / R_KJmol),
        label="MACE-MP-0",
        color="blue",
        s=10,
    )


    plt.xlabel("log10(Temperature (K))")
    plt.ylabel("log10($C_V/R$)")
    plt.title("Ice Ih heat capacity")
    plt.legend()
    # plt.savefig("heat_capacity_loglog.svg")
    plt.gca()
    return


@app.cell
def __(np):
    help(np.log)
    return


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
