import marimo

__generated_with = "0.3.3"
app = marimo.App()


@app.cell
def __():
    import marimo as mo
    return mo,


@app.cell
def __(mo):
    mo.md(
        """
        # Dispersione fononica del ghiaccio pesante

        Si simula la configurazione indicata dall'articolo Strässle2004.
        """
    )
    return


@app.cell
def __():
    import os
    import phonopy
    from ase.io import read
    from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

    small_pressure = {}

    small_pressure["basepath"] = (
        "/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/04_crystal_phonons/phonopy/MACE-ICE13-1/D2O-Ih/pressure=0.05GPa_supercell=3"
    )

    print(
        "Pressure = 0.05 GPa:",
        os.listdir(small_pressure["basepath"]),
    )

    big_pressure = {}

    big_pressure["basepath"] = (
        "/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/04_crystal_phonons/phonopy/MACE-ICE13-1/D2O-Ih/pressure=0.5GPa_supercell=3"
    )

    print(
        "Pressure = 0.5 GPa:",
        os.listdir(big_pressure["basepath"]),
    )

    for p in [small_pressure, big_pressure]:
        p["atoms"] = read(p["basepath"] + "/optimized.xyz")
        p["bandpath"] = p["atoms"].cell.bandpath()

    print("Masses:")
    print(small_pressure["atoms"].get_masses())

    labels = ["A", "G", "K", "M", "G"]
    for _p in [small_pressure, big_pressure]:
        _p["path"] = [_p["bandpath"].special_points.get(key) for key in labels]

    for _p in [small_pressure, big_pressure]:
        _p["qpoints"], _p["connections"] = get_band_qpoints_and_path_connections(
            [_p["path"]], npoints=101
        )

    for _p in [small_pressure, big_pressure]:
        _p["phonon"] = phonopy.load(_p["basepath"] + "/phonopy_params.yaml")
    return (
        big_pressure,
        get_band_qpoints_and_path_connections,
        labels,
        os,
        p,
        phonopy,
        read,
        small_pressure,
    )


@app.cell
def __(small_pressure):
    small_pressure["bandpath"].plot()
    return


@app.cell
def __(small_pressure):
    small_pressure["path"]
    return


@app.cell
def __(big_pressure):
    print(big_pressure["path"])
    return


@app.cell
def __(small_pressure):
    print(small_pressure["qpoints"])
    return


@app.cell
def __(small_pressure):
    print(small_pressure["connections"])
    return


@app.cell
def __(big_pressure, labels, small_pressure):
    for _p in [small_pressure, big_pressure]:
        _p["phonon"].run_band_structure(
            _p["qpoints"],
            path_connections=_p["connections"],
            labels=labels,
            is_legacy_plot=False,
        )
    return


@app.cell
def __():
    import matplotlib.pyplot as plt
    # from phonopy.phonon.band_structure import band_plot
    from phonopy.phonon.band_structure import BandPlot
    def band_plot(axs, frequencies, distances, path_connections, labels, fmt="r-"):
        """Return band structure plot."""
        bp = BandPlot(axs)
        bp.decorate(labels, path_connections, frequencies, distances)
        bp.plot(distances, frequencies, path_connections, fmt=fmt)
    return BandPlot, band_plot, plt


@app.cell
def __(band_plot, big_pressure, labels, plt, small_pressure):
    n = len(
        [
            x
            for x in small_pressure["phonon"]._band_structure.path_connections
            if not x
        ]
    )

    _fig, _axs = plt.subplots(
        ncols=2, nrows=1, layout="constrained", width_ratios=[9999, 1]
    )


    def _plot(phonon, connections, fmt="r-"):
        _dict = phonon.get_band_structure_dict()
        _frequencies = _dict["frequencies"]
        _distances = _dict["distances"]
        band_plot(_axs, _frequencies, _distances, connections, labels, fmt=fmt)


    _plot(small_pressure["phonon"], small_pressure["connections"], fmt="k-")
    _plot(big_pressure["phonon"], big_pressure["connections"], fmt="k--")

    # Create placeholders for the legend
    _axs[0].plot([], [], "k-", label="Pressure = 0.05 GPa")
    _axs[0].plot([], [], "k--", label="Pressure = 0.5 GPa")
    _axs[0].legend()

    _fig.delaxes(_axs[1])
    _ax = _axs[0]
    _ax.set_ylim(0, 5)
    _ax.set_ylabel("Frequency (THz)")
    _ax.grid()
    _ax.set_title("$\mathrm{D_2O}$ Ih")

    # _fig.savefig("band_structure.svg")

    _fig
    return n,


@app.cell
def __(big_pressure, plt, small_pressure):
    _fig, _ax = plt.subplots(layout="constrained")

    for _p, _fmt in zip([small_pressure, big_pressure], ["k-", "k--"]):
        # Calcola la DOS, necessario solo la prima volta, poi si può commentare
        # _p["phonon"].run_mesh([16, 16, 16])
        _p["phonon"].run_total_dos(sigma=0.05, freq_max=11)
        _p["total_dos"] = _p["phonon"].get_total_dos_dict()

        _ax.plot(_p["total_dos"]["frequency_points"], _p["total_dos"]["total_dos"], _fmt)


    def THz2K(THz):
        return THz * 47.9924307337


    def K2THz(K):
        return K * 0.0208366191233


    secax = _ax.secondary_xaxis("top", functions=(THz2K, K2THz))
    secax.set_xlabel("E (K)")

    _ax.set_ylim(bottom=0)
    _ax.set_xlim(left=0)
    _ax.set_xlabel("Frequency (THz)")
    _ax.set_ylabel("DOS (1/THz)")
    _ax.legend(["Pressure = 0.05 GPa", "Pressure = 0.5 GPa"])
    _ax.set_title("$\mathrm{D_2O}$ Ih")

    # plt.savefig("total_dos.svg")

    _fig
    return K2THz, THz2K, secax


if __name__ == "__main__":
    app.run()
