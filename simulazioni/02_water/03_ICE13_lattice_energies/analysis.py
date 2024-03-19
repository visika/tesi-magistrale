import marimo

__generated_with = "0.3.3"
app = marimo.App()


@app.cell
def __(mo):
    mo.md(
        """
        # Grafici per il confronto di MACE-MP-0 e MACE-ICE13-1 con i rispettivi pseudopotenziali di riferimento

        I dati per i 13 cristalli di ghiaccio analizzati con MACE-MP-0 medium con dispersione sono presenti nella cartella `09_ICE13_MACE`.

        Nel file `crystal_energies.csv` sono contenute le energie dei cristalli, calcolate come energia totale del sistema diviso il numero di molecole presenti, \(E_\mathrm{crystal} \coloneqq E / N_\mathrm{H_2O}\).

        Per calcolare l'energia del reticolo, si calcola la differenza \(E_\mathrm{lattice} \coloneqq E_\mathrm{crystal} - E_\mathrm{gas}\).

        La quantità \(E_\mathrm{gas}\) è stata calcolata nelle run con codice `01`. In particolare, la simulazione con migliore convergenza è nella cartella `01.8_molecule_converge_fmax_parallel_dispersion`.
        """
    )
    return


@app.cell
def __():
    import marimo as mo
    import polars as pl
    import matplotlib.pyplot as plt
    return mo, pl, plt


@app.cell
def __():
    # Leggi il valore dell'energia della molecola nello stato gassoso
    _f = open(
        "../01_molecule/01.8_molecule_converge_fmax_parallel_dispersion/medium/1e-8/e_gas.txt",
        "r",
    )
    e_gas_medium_dispersion_ev = float(_f.readlines()[0])

    _f = open(
        "../01_molecule/01.8_molecule_converge_fmax_parallel_dispersion/large/1e-8/e_gas.txt",
        "r",
    )
    e_gas_large_dispersion_ev = float(_f.readlines()[0])
    return e_gas_large_dispersion_ev, e_gas_medium_dispersion_ev


@app.cell
def __():
    # Valore di conversione tra Electron Volt per particle e KiloJoule per mole
    evp_to_kjmol = 96.4853321233
    return evp_to_kjmol,


@app.cell
def __(e_gas_medium_dispersion_ev, evp_to_kjmol, pl):
    df_medium_d = pl.read_csv("ICE13_MACE_medium/crystal_energies.csv")
    df_medium_d = df_medium_d.with_columns(
        (pl.col("e_crys") - e_gas_medium_dispersion_ev).alias("e_lattice_evp")
    )
    df_medium_d = df_medium_d.with_columns(
        (pl.col("e_lattice_evp") * evp_to_kjmol).alias("e_lattice_kjmol")
    )
    df_medium_d = df_medium_d.with_columns(
        (pl.col("e_lattice_evp") * 1e3).alias("e_lattice_mevp")
    )
    # df_medium_d
    return df_medium_d,


@app.cell
def __(e_gas_large_dispersion_ev, evp_to_kjmol, pl):
    df_large_d = pl.read_csv("ICE13_MACE_large/crystal_energies.csv")
    df_large_d = df_large_d.with_columns(
        (pl.col("e_crys") - e_gas_large_dispersion_ev).alias("e_lattice_ev")
    )
    df_large_d = df_large_d.with_columns(
        (pl.col("e_lattice_ev") * evp_to_kjmol).alias("e_lattice_kjmol")
    )
    return df_large_d,


@app.cell
def __(evp_to_kjmol, pl):
    with open(
        "/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/01_molecule/MACE-ICE13/dispersion=True/1e-8/e_gas.txt",
        "r",
    ) as _f:
        e_gas_mace_ice13_0 = float(_f.readlines()[0])

    df_mace_ice13_0 = (
        pl.read_csv("MACE-ICE13-0/crystal_energies.csv")
        .with_columns(
            (pl.col("e_crys") - e_gas_mace_ice13_0).alias("e_lattice_ev")
        )
        .with_columns(
            (pl.col("e_lattice_ev") * evp_to_kjmol).alias("e_lattice_kjmol")
        )
    )
    # df_mace_ice13_0
    return df_mace_ice13_0, e_gas_mace_ice13_0


@app.cell
def __(mo):
    mo.md("Di seguito la struttura dati dei polimorfi del ghiaccio analizzati con MACE-ICE13-1, modello addestrato su revPBE-D3.")
    return


@app.cell(hide_code=True)
def __(evp_to_kjmol, pl):
    with open(
        "/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/01_molecule/MACE-ICE13-1/e_gas.txt",
        "r",
    ) as _f:
        e_gas_mace_ice13_1 = float(_f.readlines()[0])

    df_mace_ice13_1 = pl.read_csv("MACE-ICE13-1/crystal_energies.csv")
    df_mace_ice13_1 = df_mace_ice13_1.with_columns(
        (pl.col("e_crys") - e_gas_mace_ice13_1).alias("e_lattice_ev")
    )
    df_mace_ice13_1 = df_mace_ice13_1.with_columns(
        (pl.col("e_lattice_ev") * evp_to_kjmol).alias("e_lattice_kjmol")
    )

    with open(
        "/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/01_molecule/MACE-ICE13-1-D3/e_gas.txt",
        "r",
    ) as _f:
        e_gas_mace_ice13_1_d3 = float(_f.readlines()[0])

    df_mace_ice13_1
    return df_mace_ice13_1, e_gas_mace_ice13_1, e_gas_mace_ice13_1_d3


@app.cell(hide_code=True)
def __(evp_to_kjmol, pl):
    dmc_mev_table = [
        ["Ih", -616],
        ["II", -613],
        ["III", -603],
        ["IV", -576],
        ["VI", -597],
        ["VII", -564],
        ["VIII", -582],
        ["IX", -610],
        ["XI", -614],
        ["XIII", -594],
        ["XIV", -598],
        ["XV", -598],
        ["XVII", -598],
    ]

    pbe_d3_kjmol_table = [
        ["Ih", -70.80],
        ["II", -68.37],
        ["III", -68.47],
        ["IV", -65.82],
        ["VI", -66.17],
        ["VII", -60.66],
        ["VIII", -61.79],
        ["IX", -68.97],
        ["XI", -71.29],
        ["XIII", -67.76],
        ["XIV", -67.11],
        ["XV", -66.00],
        ["XVII", -69.75],
    ]

    df_dmc = pl.DataFrame(dmc_mev_table, schema=["structure", "e_lattice_mevp"])
    df_dmc = df_dmc.with_columns(
        (pl.col("e_lattice_mevp") * evp_to_kjmol / 1000).alias("e_lattice_kjmol")
    )

    df_pbe_d3 = pl.DataFrame(
        pbe_d3_kjmol_table, schema=["structure", "e_lattice_kjmol"]
    )
    return df_dmc, df_pbe_d3, dmc_mev_table, pbe_d3_kjmol_table


@app.cell(hide_code=True)
def __(pl):
    rev_pbe_d3_kjmol = [
        ["Ih", -59.01],
        ["II", -57.75],
        ["III", -56.69],
        ["IV", -55.03],
        ["VI", -56.16],
        ["VII", -54.83],
        ["VIII", -55.74],
        ["IX", -57.41],
        ["XI", -59.25],
        ["XIII", -56.71],
        ["XIV", -56.30],
        ["XV", -56.07],
        ["XVII", -58.00],
    ]

    df_rev_pbe_d3 = pl.DataFrame(
        rev_pbe_d3_kjmol, schema=["structure", "e_lattice_kjmol"]
    )
    return df_rev_pbe_d3, rev_pbe_d3_kjmol


@app.cell(hide_code=True)
def __(
    df_dmc,
    df_large_d,
    df_mace_ice13_1,
    df_medium_d,
    df_pbe_d3,
    df_rev_pbe_d3,
    evp_to_kjmol,
    plt,
):
    # Diffusion Monte Carlo
    plt.plot(
        df_dmc["structure"],
        df_dmc["e_lattice_mevp"] * evp_to_kjmol / 1000,
        marker="*",
        label="DMC",
        markersize=10,
        linestyle="-.",
    )

    # revPBE-D3
    plt.plot(
        df_rev_pbe_d3["structure"],
        df_rev_pbe_d3["e_lattice_kjmol"],
        label="revPBE-D3",
        marker="s",
        linestyle="--",
    )

    # MACE-ICE13-1
    plt.plot(
        df_mace_ice13_1["structure"],
        df_mace_ice13_1["e_lattice_kjmol"],
        marker="o",
        label="MACE-ICE13-1-D3",
    )

    # PBE-D3
    plt.plot(
        df_pbe_d3["structure"],
        df_pbe_d3["e_lattice_kjmol"],
        marker="s",
        label="PBE-D3",
        linestyle="--",
    )

    # MACE-MP-0 medium
    plt.plot(
        df_medium_d["structure"],
        df_medium_d["e_lattice_evp"] * evp_to_kjmol,
        marker="o",
        label="MACE-MP-0 medium D3",
    )

    # MACE-MP-0 large
    plt.plot(
        df_large_d["structure"],
        df_large_d["e_lattice_kjmol"],
        marker="o",
        label="MACE-MP-0 large D3",
    )

    # plt.plot(
    #     df_mace_d_opt["structure"],
    #     df_mace_d_opt["e_lattice_kjmol"],
    #     marker="o",
    #     label="MACE medium-D opt",
    # )

    # MACE-ICE13-0
    # plt.plot(
    #     df_mace_ice13_0["structure"],
    #     df_mace_ice13_0["e_lattice_kjmol"],
    #     marker="o",
    #     label="MACE-ICE13-0 D3",
    # )

    plt.xlabel("Structure")
    plt.ylabel("Lattice energy (kJ/mol)")
    plt.grid()
    plt.legend(loc="lower right")
    plt.title("Absolute lattice energy")
    return


@app.cell
def __(pl):
    from os.path import join

    def get_the_dataframe(path):
        folders = [
            "Ih",
            "II",
            "III",
            "IV",
            "VI",
            "VII",
            "VIII",
            "IX",
            "XI",
            "XIII",
            "XIV",
            "XV",
            "XVII",
        ]

        table = []

        for folder in folders:
            file_path = join(path, folder, "e_crys.txt")
            with open(file_path, "r") as file:
                value = float(file.read())
                table.append([folder, value])

        df = pl.DataFrame(table, schema=["structure", "e_crys"])
        return df
    return get_the_dataframe, join


@app.cell
def __(mo):
    mo.md("# Relative lattice energy")
    return


@app.cell(hide_code=True)
def __(
    df_dmc,
    df_mace_ice13_1,
    df_medium_d,
    df_pbe_d3,
    df_rev_pbe_d3,
    pl,
    plt,
):
    def get_relative_lattice_energy(dataframe):
        # Select lattice energy in the row with structure Ih
        baseline = dataframe.filter(pl.col("structure") == "Ih")[
            "e_lattice_kjmol"
        ][0]
        # Calculate the difference between the lattice energy of the other structures and the lattice energy of the structure Ih
        relative = dataframe.with_columns(
            (pl.col("e_lattice_kjmol") - baseline).alias(
                "e_lattice_relative_kjmol"
            )
        )
        return relative


    _dmc = get_relative_lattice_energy(df_dmc)
    _rev_pbe_d3 = get_relative_lattice_energy(df_rev_pbe_d3)
    _mace_ice13_1 = get_relative_lattice_energy(df_mace_ice13_1)
    _pbe_d3 = get_relative_lattice_energy(df_pbe_d3)
    _mace_mp_0_medium = get_relative_lattice_energy(df_medium_d)

    # Plot the relative lattice energy

    plt.plot(
        _dmc["structure"],
        _dmc["e_lattice_relative_kjmol"],
        marker="*",
        label="DMC",
        markersize=10,
        linestyle="-.",
    )

    plt.plot(
        _rev_pbe_d3["structure"],
        _rev_pbe_d3["e_lattice_relative_kjmol"],
        label="revPBE-D3",
        marker="s",
        linestyle="--",
        linewidth=3,
    )

    plt.plot(
        _mace_ice13_1["structure"],
        _mace_ice13_1["e_lattice_relative_kjmol"],
        marker="o",
        label="MACE-ICE13-1-D3",
        linewidth=1,
        markersize=4,
    )

    plt.plot(
        _pbe_d3["structure"],
        _pbe_d3["e_lattice_relative_kjmol"],
        marker="s",
        label="PBE-D3",
        linestyle="--",
    )

    plt.plot(
        _mace_mp_0_medium["structure"],
        _mace_mp_0_medium["e_lattice_relative_kjmol"],
        marker="o",
        label="MACE-MP-0 medium D3",
        linestyle="-",
    )

    plt.xlabel("Structure")
    plt.ylabel("Relative lattice energy (kJ/mol)")
    plt.grid()

    plt.legend()
    return get_relative_lattice_energy,


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        """
        # Scatter plot tra MACE-MP-0 e PBE-D3

        Lo scatter plot mostra che effettivamente la geometria non ottimizzata fornisce risultati più vicini al modello di riferimento.

        Si osserva uno scarto pressocché costante nei valori tra il modello in esame e il modello di riferimento.

        Il sospetto di Zen, e ora anche il mio, ricadono sulla quantità costante all'interno della relazione che definisce l'energia del reticolo, l'energia della molecola nello stato gassoso:

        \[ E_\mathrm{lattice} \coloneqq E_\mathrm{crystal} - E_\mathrm{gas} \]
        """
    )
    return


@app.cell(hide_code=True)
def __(df_large_d, df_medium_d, df_pbe_d3, plt):
    # Make a scatter plot between MACE medium-D and PBE-D3

    # Make the figure a square
    plt.figure(figsize=(6, 6))

    plt.scatter(
        df_medium_d["e_lattice_kjmol"],
        df_pbe_d3["e_lattice_kjmol"],
        label="MACE-MP-0 medium D3",
        zorder=2,
        # s=15,
        marker="s",
    )

    plt.scatter(
        df_large_d["e_lattice_kjmol"],
        df_pbe_d3["e_lattice_kjmol"],
        label="MACE-MP-0 large D3",
        # s=15,
        marker="o",
    )

    plt.xlabel("MACE-MP-0 D3 (kJ/mol)")
    plt.ylabel("PBE-D3 (kJ/mol)")

    plt.axline((0, 0), slope=1, linestyle="--")

    plt.xlim(-80, -60)
    plt.ylim(-80, -60)

    plt.grid()
    plt.legend()
    # plt.savefig("scatterplot_mace_vs_pbe.png", dpi=300)
    return


@app.cell
def __(mo):
    mo.md("# Scatter plot tra MACE-ICE13-1 e revPBE-D3")
    return


@app.cell(hide_code=True)
def __(df_mace_ice13_1, df_rev_pbe_d3, plt):
    plt.figure(figsize=(6, 6))

    plt.scatter(
        df_mace_ice13_1["e_lattice_kjmol"],
        df_rev_pbe_d3["e_lattice_kjmol"],
        label="MACE-ICE13-1 D3",
    )

    plt.xlabel("MACE-ICE13-1 D3 (kJ/mol)")
    plt.ylabel("revPBE-D3 (kJ/mol)")

    plt.axline((0, 0), slope=1, linestyle="--")

    plt.xlim(-60, -54)
    plt.ylim(-60, -54)

    plt.grid()
    plt.legend()
    return


if __name__ == "__main__":
    app.run()
