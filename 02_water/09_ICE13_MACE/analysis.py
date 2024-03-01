import marimo

__generated_with = "0.2.13"
app = marimo.App(width="full")


@app.cell
def __():
    import marimo as mo
    import polars as pl
    return mo, pl


@app.cell
def __(e_gas, evp_to_kjmol, pl):
    df = pl.read_csv("crystal_energies.csv")
    df = df.with_columns((pl.col("e_crys") - e_gas).alias("e_lattice"))
    df = df.with_columns((pl.col("e_lattice") * 1e3).alias("e_lattice_mev"))
    df = df.with_columns(
        (pl.col("e_lattice") * evp_to_kjmol).alias("e_lattice_kjmol")
    )
    df
    return df,


app._unparsable_cell(
    r"""
    df_large_d = pl.read_csv(\"../09.1_ICE13_MACE_large/crystal_energies.csv\")
    df_large_d = df_large_d.with_columns((pl.col(\"e_crys\") - e_gas_large_d).alias(\"e_lattice\")))
    """,
    name="__"
)


@app.cell
def __(e_gas_large_d):
    e_gas = -14.170072284496875
    e_gas_large_d
    return e_gas,


@app.cell
def __():
    evp_to_kjmol = 96.4916
    return evp_to_kjmol,


@app.cell
def __():
    dmc_table = [
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

    pbe_d3_kj_table = [
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
    return dmc_table, pbe_d3_kj_table


@app.cell
def __(dmc_table, pbe_d3_kj_table, pl):
    df_dmc = pl.DataFrame(dmc_table, schema=["structure", "e_lattice"])
    df_pbe_d3 = pl.DataFrame(
        pbe_d3_kj_table, schema=["structure", "e_lattice_kjmol"]
    )
    return df_dmc, df_pbe_d3


@app.cell
def __():
    return


@app.cell
def __(df, df_dmc, df_mace_d_opt, df_pbe_d3):
    df_join = df.join(
        df_dmc, left_on="structure", right_on="structure", suffix="_dmc"
    )
    df_join = df_join.join(
        df_pbe_d3, left_on="structure", right_on="structure", suffix="_pbe_d3"
    )
    df_join = df_join.join(
        df_mace_d_opt, left_on="structure", right_on="structure", suffix="_mace_opt"
    )
    df_join
    return df_join,


@app.cell
def __(df_join):
    df_join.plot("e_lattice_mev", "e_lattice_dmc", kind="scatter")
    return


@app.cell
def __():
    import matplotlib.pyplot as plt
    return plt,


@app.cell
def __(df_join, df_mace_d_opt, df_pbe_d3, evp_to_kjmol, plt):
    plt.plot(
        df_join["structure"],
        df_join["e_lattice_mev"] * evp_to_kjmol / 1000,
        marker="o",
        label="MACE medium-D no opt",
    )

    plt.plot(
        df_join["structure"],
        df_join["e_lattice_dmc"] * evp_to_kjmol / 1000,
        marker="*",
        label="DMC",
    )

    plt.plot(
        df_pbe_d3["structure"],
        df_pbe_d3["e_lattice_kjmol"],
        marker="^",
        label="PBE-D3",
    )

    plt.plot(
        df_mace_d_opt["structure"],
        df_mace_d_opt["e_lattice_kjmol"],
        marker="o",
        label="MACE medium-D opt",
    )

    plt.xlabel("Structure")
    plt.ylabel("Lattice energy (kJ/mol)")
    plt.grid()
    plt.legend()
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
def __(e_gas, evp_to_kjmol, get_the_dataframe, pl):
    df_mace_d_opt = get_the_dataframe("../06_crystal_structures_mace_dispersion")
    df_mace_d_opt = df_mace_d_opt.with_columns(
        ((pl.col("e_crys") - e_gas) * evp_to_kjmol).alias("e_lattice_kjmol")
    )
    return df_mace_d_opt,


@app.cell
def __(mo):
    mo.md(
        """
        # Scatter plot tra MACE e PBE-D3

        Lo scatter plot mostra che effettivamente la geometria non ottimizzata fornisce risultati più vicini al modello di riferimento.

        Si osserva uno scarto pressocché costante nei valori tra il modello in esame e il modello di riferimento.

        Il sospetto di Zen, e ora anche il mio, ricadono sulla quantità costante all'interno della relazione che definisce l'energia del reticolo, l'energia della molecola nello stato gassoso:

        \[ E_\mathrm{lattice} \coloneqq E_\mathrm{crystal} - E_\mathrm{gas} \]
        """
    )
    return


@app.cell
def __(df_join, plt):
    # Make a scatter plot between MACE medium-D no opt and PBE-D3
    plt.scatter(
        df_join["e_lattice_kjmol"],
        df_join["e_lattice_kjmol_pbe_d3"],
        label="no opt",
        zorder=2,
    )
    plt.scatter(
        df_join["e_lattice_kjmol_mace_opt"],
        df_join["e_lattice_kjmol_pbe_d3"],
        label="opt",
        s=15,
    )
    plt.xlabel("MACE-MP-0 medium-D (kJ/mol)")
    plt.ylabel("PBE-D3 (kJ/mol)")

    plt.axline((0, 0), slope=1, linestyle="--")

    plt.xlim(-80, -60)
    plt.ylim(-80, -60)

    plt.grid()
    plt.legend()
    # plt.savefig("scatterplot_mace_vs_pbe.png", dpi=300)
    return


if __name__ == "__main__":
    app.run()
