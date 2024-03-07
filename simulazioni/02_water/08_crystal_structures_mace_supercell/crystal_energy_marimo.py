import marimo

__generated_with = "0.2.13"
app = marimo.App(width="full")


@app.cell
def __(merged):
    # Define the sorting order
    sorting_order = [
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
    # Sort the dataframe
    sorted = merged.set_index("folder").loc[sorting_order].reset_index()
    return sorted, sorting_order


@app.cell
def __(plt, sorted):
    plt.figure(figsize=(10, 10))
    plt.plot(sorted["e_lattice_dmc_kj"], label="DMC", marker="*", linestyle="-")

    plt.plot(
        sorted["e_lattice_medium_kj"],
        label="MACE medium dispersion",
        marker="s",
        linestyle="-",
    )

    plt.plot(
        sorted["e_lattice_without_dispersion_kj"],
        label="MACE medium no dispersion",
        marker="^",
        linestyle="-",
    )

    plt.plot(
        sorted["e_lattice_b3lypd4_kj"], label="B3LYP-D4", marker="o", linestyle="-"
    )

    plt.plot(sorted["e_lattice_n2p2_kj"], label="n2p2", marker="x", linestyle="-")
    plt.plot(sorted["e_lattice_lda_kj"], label="LDA", marker="v", linestyle="-")
    plt.plot(
        sorted["e_lattice_only_geometry_kj"],
        label="MACE medium dispersion\nonly geometry optimization",
        marker="d",
        linestyle="-",
    )

    plt.plot(
        sorted["e_lattice_pbe_d3_kj"],
        label="PBE-D3",
        marker="p",
        linestyle="-",
    )

    plt.plot(
        sorted["e_lattice_rev_pbe0_d3_kj"],
        label="revPBE0-D3",
        marker="h",
        linestyle="-",
    )

    plt.plot(
        sorted["e_lattice_supercell_kj"],
        label="MACE medium supercell 2x2x2",
        marker="8",
        linestyle="-",
    )

    plt.plot(
        sorted["e_lattice_no_dispersion_large_kj"],
        label="MACE no dispersion large",
        marker="P",
        linestyle="-",
    )

    plt.xlabel("Crystal structure")
    plt.ylabel("Absolute lattice energy (kJ/mol)")
    plt.xticks(range(len(sorted["folder"])), sorted["folder"])

    # Put the legend outside
    plt.legend(loc="lower center", bbox_to_anchor=(0.5, -0.2), ncol=3)

    plt.grid()
    # Make all linewidths smaller
    [line.set_linewidth(0.5) for line in plt.gca().lines]
    # Make all symbols smaller
    [line.set_markersize(5) for line in plt.gca().lines]
    plt.title("Absolute lattice energy")

    plt.tight_layout()
    # plt.savefig("absolute_lattice_energy.png", dpi=300)
    plt.show()
    return


@app.cell(hide_code=True)
def __(plt, sorted):
    plt.figure(figsize=(10, 10))
    plt.plot(
        sorted["e_lattice_dmc_kj_rel"], label="DMC", marker="*", linestyle="-"
    )

    plt.plot(
        sorted["e_lattice_medium_kj_rel"],
        label="MACE medium dispersion",
        marker="s",
        linestyle="-",
    )

    plt.plot(
        sorted["e_lattice_without_dispersion_kj_rel"],
        label="MACE medium no dispersion",
        marker="^",
        linestyle="-",
    )

    plt.plot(
        sorted["e_lattice_b3lypd4_kj_rel"],
        label="B3LYP-D4",
        marker="o",
        linestyle="-",
    )

    plt.plot(
        sorted["e_lattice_n2p2_kj_rel"], label="n2p2", marker="x", linestyle="-"
    )

    plt.plot(
        sorted["e_lattice_lda_kj_rel"], label="LDA", marker="v", linestyle="-"
    )

    plt.plot(
        sorted["e_lattice_only_geometry_kj_rel"],
        label="MACE medium dispersion\nonly geometry optimization",
        marker="d",
        linestyle="-",
    )

    plt.plot(
        sorted["e_lattice_pbe_d3_kj_rel"],
        label="PBE-D3",
        marker="p",
        linestyle="-",
    )

    plt.plot(
        sorted["e_lattice_rev_pbe0_d3_kj_rel"],
        label="revPBE0-D3",
        marker="h",
        linestyle="-",
    )

    plt.plot(
        sorted["e_lattice_supercell_kj_rel"],
        label="MACE medium supercell 2x2x2",
        marker="8",
        linestyle="-",
    )

    plt.plot(
        sorted["e_lattice_no_dispersion_large_kj_rel"],
        label="MACE no dispersion large",
        marker="P",
        linestyle="-",
    )

    plt.xlabel("Crystal structure")
    plt.ylabel("Relative lattice energy (kJ/mol)")
    plt.xticks(range(len(sorted["folder"])), sorted["folder"])

    # Put the legend outside
    plt.legend(loc="lower center", bbox_to_anchor=(0.5, -0.2), ncol=3)

    plt.grid()
    # Make all linewidths smaller
    [line.set_linewidth(0.5) for line in plt.gca().lines]
    # Make all symbols smaller
    [line.set_markersize(5) for line in plt.gca().lines]
    plt.title("Relative lattice energy")

    plt.tight_layout()
    # plt.savefig("relative_lattice_energy.png", dpi=300)
    plt.show()
    return


@app.cell
def __():
    # Energia della molecola con dispersione in eV
    e_gas = -14.170072284496875
    # Energia della molecola senza dispersione in eV
    e_gas_without_dispersion = -14.160290476990022
    return e_gas, e_gas_without_dispersion


@app.cell(hide_code=True)
def __():
    # Import listdir
    import os
    from os import listdir
    from os.path import isfile, join, isdir

    import pandas as pd
    import matplotlib.pyplot as plt

    os.getcwd()
    return isdir, isfile, join, listdir, os, pd, plt


@app.cell
def __(isdir, join, listdir, pd):
    def get_the_dataframe(path):
        # List folders
        folders = [f for f in listdir(path) if isdir(join(path, f))]
        # List folders with roman numerals
        folders_roman = [f for f in folders if f[0].isupper()]

        # Put the value contained in the file e_crys.txt from each directory in a table
        # Define the table
        table = []
        # Loop over the folders
        for folder in folders_roman:
            # Define the file path
            file_path = join(path, folder, "e_crys.txt")
            # Open the file
            with open(file_path, "r") as file:
                # Read the value
                value = float(file.read())
                # Append the value to the table
                table.append([folder, value])

        # Render the table as a pandas dataframe
        df = pd.DataFrame(table, columns=["folder", "e_crys"])
        return df
    return get_the_dataframe,


@app.cell
def __(e_gas, get_the_dataframe):
    df_supercell = get_the_dataframe(".")
    df_supercell["e_lattice_supercell"] = df_supercell["e_crys"] - e_gas
    df_supercell
    return df_supercell,


@app.cell
def __(e_gas, get_the_dataframe):
    df = get_the_dataframe("../06_crystal_structures_mace_dispersion")
    # Calculate the difference of the values in the column e_crys with respect to the value associated to the folder Ih
    # Define the reference value
    ref_value = df[df["folder"] == "Ih"]["e_crys"].values[0]
    # Calculate the lattice energy as the difference between e_crys and e_gas
    df["e_lattice"] = df["e_crys"] - e_gas
    df
    return df, ref_value


@app.cell
def __(e_gas_without_dispersion, get_the_dataframe):
    df_without_dispersion = get_the_dataframe(
        "../04_crystal_structures_mace/medium"
    )
    df_without_dispersion["e_lattice_without_dispersion"] = (
        df_without_dispersion["e_crys"] - e_gas_without_dispersion
    )
    return df_without_dispersion,


@app.cell
def __(e_gas_without_dispersion, get_the_dataframe):
    df_no_dispersion_large = get_the_dataframe(
        "../04_crystal_structures_mace/large"
    )
    df_no_dispersion_large["e_lattice_no_dispersion_large"] = (
        df_no_dispersion_large["e_crys"] - e_gas_without_dispersion
    )
    return df_no_dispersion_large,


@app.cell
def __(get_the_dataframe):
    # Energia con n2p2
    e_gas_n2p2 = -468.46641416393885
    df_n2p2 = get_the_dataframe("../03_crystal_structures_n2p2")
    df_n2p2["e_lattice_n2p2"] = df_n2p2["e_crys"] - e_gas_n2p2
    return df_n2p2, e_gas_n2p2


@app.cell
def __(e_gas, get_the_dataframe):
    df_only_geometry = get_the_dataframe(
        "../07_crystal_structures_mace_dispersion_only_geometry_optimization"
    )
    df_only_geometry["e_lattice_only_geometry"] = (
        df_only_geometry["e_crys"] - e_gas
    )
    return df_only_geometry,


@app.cell(hide_code=True)
def __(pd):
    # Build the DMC table by hand. Define a dataframe with columns 'folder', 'e_lattice'.
    # Define the table
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
    df_dmc = pd.DataFrame(dmc_table, columns=["folder", "e_lattice"])
    return df_dmc, dmc_table


@app.cell(hide_code=True)
def __(pd):
    b3lypd4_table = [
        ["Ih", -653],
        ["II", -633],
        ["III", -626],
        ["IV", -602],
        ["VI", -608],
        ["VII", -557],
        ["VIII", -570],
        ["IX", -636],
        ["XI", -655],
        ["XIII", -624],
        ["XIV", -618],
        ["XV", -608],
        ["XVII", -643],
    ]
    df_b3lypd4 = pd.DataFrame(b3lypd4_table, columns=["folder", "e_lattice_b3lypd4"])
    return b3lypd4_table, df_b3lypd4


@app.cell(hide_code=True)
def __(pd):
    lda_table = [
        ["Ih", -1037],
        ["II", -978],
        ["III", -994],
        ["IV", -943],
        ["VI", -943],
        ["VII", -868],
        ["VIII", -877],
        ["IX", -988],
        ["XI", -1049],
        ["XIII", -964],
        ["XIV", -953],
        ["XV", -936],
        ["XVII", -1033],
    ]
    df_lda = pd.DataFrame(lda_table, columns=["folder", "e_lattice_lda"])
    return df_lda, lda_table


@app.cell
def __(pd):
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
    df_pbe_d3_kj = pd.DataFrame(
        pbe_d3_kj_table, columns=["folder", "e_lattice_pbe_d3_kj"]
    )
    return df_pbe_d3_kj, pbe_d3_kj_table


@app.cell(hide_code=True)
def __(pd):
    rev_pbe0_d3_kj_table = [
        ["Ih", -58.45],
        ["II", -57.58],
        ["III", -55.35],
        ["IV", -54.47],
        ["VI", -56.20],
        ["VII", -55.61],
        ["VIII", -56.67],
        ["IX", -56.61],
        ["XI", -57.51],
        ["XIII", -56.33],
        ["XIV", -56.08],
        ["XV", -56.15],
        ["XVII", -56.44],
    ]
    df_rev_pbe0_d3_kj = pd.DataFrame(
        rev_pbe0_d3_kj_table, columns=["folder", "e_lattice_rev_pbe0_d3_kj"]
    )
    return df_rev_pbe0_d3_kj, rev_pbe0_d3_kj_table


@app.cell
def __():
    # Calcola il fattore di conversione tra kJ/mol ed eV
    # DMC in tabella II per il lattice energy di Ih: -59.45(7) kJ/mol
    dmc_kj = -59.45
    dmc_ev = -616
    conversion_factor = dmc_kj / dmc_ev
    conversion_factor
    return conversion_factor, dmc_ev, dmc_kj


@app.cell
def __(
    conversion_factor,
    df,
    df_b3lypd4,
    df_dmc,
    df_lda,
    df_n2p2,
    df_no_dispersion_large,
    df_only_geometry,
    df_pbe_d3_kj,
    df_rev_pbe0_d3_kj,
    df_supercell,
    df_without_dispersion,
    pd,
):
    merged = pd.merge(df, df_dmc, on="folder", suffixes=("_medium", "_dmc"))
    merged = pd.merge(merged, df_b3lypd4, on="folder")
    merged = pd.merge(
        merged,
        df_without_dispersion,
        on="folder",
        suffixes=("_merged", "_without_dispersion"),
    )
    merged = pd.merge(merged, df_n2p2, on="folder", suffixes=("_merged", "_n2p2"))
    merged = pd.merge(merged, df_lda, on="folder", suffixes=("_merged", "_lda"))
    merged = pd.merge(merged, df_only_geometry, on="folder")

    merged["e_lattice_medium"] = merged["e_lattice_medium"] * 1e3
    merged["e_lattice_without_dispersion"] = (
        merged["e_lattice_without_dispersion"] * 1e3
    )

    merged = pd.merge(
        merged,
        df_no_dispersion_large,
        on="folder",
        suffixes=("_merged", "_no_dispersion_large"),
    )
    merged["e_lattice_no_dispersion_large"] = (
        merged["e_lattice_no_dispersion_large"] * 1e3
    )
    # Convert e_lattice_no_dispersion_large to kJ/mol
    merged["e_lattice_no_dispersion_large_kj"] = (
        merged["e_lattice_no_dispersion_large"] * conversion_factor
    )
    # Compute the relative lattice energy
    merged["e_lattice_no_dispersion_large_kj_rel"] = (
        merged["e_lattice_no_dispersion_large_kj"]
        - merged[merged["folder"] == "Ih"][
            "e_lattice_no_dispersion_large_kj"
        ].values[0]
    )

    merged["e_lattice_n2p2"] = merged["e_lattice_n2p2"] * 1e3
    merged["e_lattice_only_geometry"] = merged["e_lattice_only_geometry"] * 1e3
    # Calcola il lattice energy in kJ/mol
    merged["e_lattice_dmc_kj"] = merged["e_lattice_dmc"] * conversion_factor
    # Select the row with Ih
    dmc_Ih = merged[merged["folder"] == "Ih"]["e_lattice_dmc_kj"].values[0]
    # Compute the relative lattice energy
    merged["e_lattice_dmc_kj_rel"] = merged["e_lattice_dmc_kj"] - dmc_Ih

    # Convert e_lattice_medium to kJ/mol
    merged["e_lattice_medium_kj"] = merged["e_lattice_medium"] * conversion_factor
    # Select the row with Ih
    medium_Ih = merged[merged["folder"] == "Ih"]["e_lattice_medium_kj"].values[0]
    # Compute the relative lattice energy for medium
    merged["e_lattice_medium_kj_rel"] = merged["e_lattice_medium_kj"] - medium_Ih

    merged["e_lattice_b3lypd4_kj"] = (
        merged["e_lattice_b3lypd4"] * conversion_factor
    )
    # Relative
    merged["e_lattice_b3lypd4_kj_rel"] = (
        merged["e_lattice_b3lypd4_kj"]
        - merged[merged["folder"] == "Ih"]["e_lattice_b3lypd4_kj"].values[0]
    )

    # Convert e_lattice_lda to kJ/mol
    merged["e_lattice_lda_kj"] = merged["e_lattice_lda"] * conversion_factor
    # Relative
    merged["e_lattice_lda_kj_rel"] = (
        merged["e_lattice_lda_kj"]
        - merged[merged["folder"] == "Ih"]["e_lattice_lda_kj"].values[0]
    )

    # Convert e_lattice_without_dispersion to kJ/mol
    merged["e_lattice_without_dispersion_kj"] = (
        merged["e_lattice_without_dispersion"] * conversion_factor
    )
    # Select the row with Ih
    without_dispersion_Ih = merged[merged["folder"] == "Ih"][
        "e_lattice_without_dispersion_kj"
    ].values[0]
    # Compute the relative lattice energy for without dispersion
    merged["e_lattice_without_dispersion_kj_rel"] = (
        merged["e_lattice_without_dispersion_kj"] - without_dispersion_Ih
    )

    # Convert n2p2 to kJ/mol
    merged["e_lattice_n2p2_kj"] = merged["e_lattice_n2p2"] * conversion_factor
    # Select the row with Ih
    n2p2_Ih = merged[merged["folder"] == "Ih"]["e_lattice_n2p2_kj"].values[0]
    # Compute the relative lattice energy for n2p2
    merged["e_lattice_n2p2_kj_rel"] = merged["e_lattice_n2p2_kj"] - n2p2_Ih

    # Convert e_lattice_only_geometry to kJ/mol
    merged["e_lattice_only_geometry_kj"] = (
        merged["e_lattice_only_geometry"] * conversion_factor
    )

    # Compute the relative lattice energy for only geometry
    merged["e_lattice_only_geometry_kj_rel"] = (
        merged["e_lattice_only_geometry_kj"]
        - merged[merged["folder"] == "Ih"]["e_lattice_only_geometry_kj"].values[0]
    )

    ### PBE-D3
    merged = pd.merge(merged, df_pbe_d3_kj, on="folder")
    # Compute the relative lattice energy for PBE-D3
    merged["e_lattice_pbe_d3_kj_rel"] = (
        merged["e_lattice_pbe_d3_kj"]
        - merged[merged["folder"] == "Ih"]["e_lattice_pbe_d3_kj"].values[0]
    )

    merged = pd.merge(merged, df_rev_pbe0_d3_kj, on="folder")
    # Compute the relative lattice energy for revPBE0-D3
    merged["e_lattice_rev_pbe0_d3_kj_rel"] = (
        merged["e_lattice_rev_pbe0_d3_kj"]
        - merged[merged["folder"] == "Ih"]["e_lattice_rev_pbe0_d3_kj"].values[0]
    )

    ### Supercell 2x2x2
    merged = pd.merge(
        merged, df_supercell, on="folder", suffixes=(None, "_supercell")
    )
    # Convert eV to meV
    merged["e_lattice_supercell"] = merged["e_lattice_supercell"] * 1e3
    # Convert e_lattice_supercell to kJ/mol
    merged["e_lattice_supercell_kj"] = (
        merged["e_lattice_supercell"] * conversion_factor
    )
    # Compute the relative lattice energy for supercell
    merged["e_lattice_supercell_kj_rel"] = (
        merged["e_lattice_supercell_kj"]
        - merged[merged["folder"] == "Ih"]["e_lattice_supercell_kj"].values[0]
    )

    merged
    return dmc_Ih, medium_Ih, merged, n2p2_Ih, without_dispersion_Ih


@app.cell
def __(df_supercell):
    df_supercell
    return


@app.cell
def __(merged):
    # Compute the mean absolute error
    mae_medium = (merged["e_lattice_medium"] - merged["e_lattice_dmc"]).abs().mean()
    mae_without_dispersion = (
        (merged["e_lattice_without_dispersion"] - merged["e_lattice_dmc"]).abs().mean()
    )
    mae_b3lypd4 = (merged["e_lattice_b3lypd4"] - merged["e_lattice_dmc"]).abs().mean()
    mae_n2p2 = (merged["e_lattice_n2p2"] - merged["e_lattice_dmc"]).abs().mean()
    mae_medium, mae_without_dispersion, mae_b3lypd4, mae_n2p2
    return mae_b3lypd4, mae_medium, mae_n2p2, mae_without_dispersion


@app.cell
def __(merged):
    mae_pbe_d3_kj = (merged["e_lattice_pbe_d3_kj"] - merged["e_lattice_dmc_kj"]).abs().mean()
    mae_pbe_d3_kj
    return mae_pbe_d3_kj,


@app.cell
def __(mae_b3lypd4, mae_medium, mae_n2p2, mae_without_dispersion):
    def get_variable_name(variable):
        variable_name = [
            name for name, value in globals().items() if value is variable
        ][0]
        return variable_name

    for e in [mae_medium, mae_without_dispersion, mae_b3lypd4, mae_n2p2]:
        print(f"{get_variable_name(e)}: {e}")
    return e, get_variable_name


@app.cell
def __(merged):
    merged.columns
    return


@app.cell
def __(merged):
    (merged["e_lattice_medium_kj"] - merged["e_lattice_supercell_kj"]).mean()
    return


if __name__ == "__main__":
    app.run()
