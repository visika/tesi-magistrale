import marimo

__generated_with = "0.2.8"
app = marimo.App()


@app.cell
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
    # Energia della molecola con dispersione
    e_gas = -14.170072284496875
    # Carica i valori del modello medium con dispersione
    data_path = "."
    # List folders
    folders = [f for f in listdir(data_path) if isdir(join(data_path, f))]
    # List folders with roman numerals
    folders_roman = [f for f in folders if f[0].isupper()]

    # Put the value contained in the file e_crys.txt from each directory in a table
    # Define the table
    table = []
    # Loop over the folders
    for folder in folders_roman:
        # Define the file path
        file_path = join(data_path, folder, "e_crys.txt")
        # Open the file
        with open(file_path, "r") as file:
            # Read the value
            value = float(file.read())
            # Append the value to the table
            table.append([folder, value])

    # Render the table as a pandas dataframe
    df = pd.DataFrame(table, columns=["folder", "e_crys"])

    # Calculate the difference of the values in the column e_crys with respect to the value associated to the folder Ih
    # Define the reference value
    ref_value = df[df["folder"] == "Ih"]["e_crys"].values[0]
    # Calculate the difference
    df["diff"] = df["e_crys"] - ref_value

    # Calculate the lattice energy as the difference between e_crys and e_gas
    df["e_lattice"] = df["e_crys"] - e_gas
    df
    return (
        data_path,
        df,
        e_gas,
        file,
        file_path,
        folder,
        folders,
        folders_roman,
        ref_value,
        table,
        value,
    )


@app.cell
def __(isdir, join, listdir, pd):
    # Energia della molecola senza dispersione
    e_gas_without_dispersion = -14.160290476990022
    # Carica i valori del modello medium senza dispersione
    data_path_without_dispersion = "../04_crystal_structures_mace/medium"
    # List folders
    folders_without_dispersion = [
        f
        for f in listdir(data_path_without_dispersion)
        if isdir(join(data_path_without_dispersion, f))
    ]
    # List folders with roman numerals
    folders_roman_without_dispersion = [
        f for f in folders_without_dispersion if f[0].isupper()
    ]
    # Put the value contained in the file e_crys.txt from each directory in a table
    # Define the table
    table_without_dispersion = []
    # Loop over the folders
    for folder_without_dispersion in folders_roman_without_dispersion:
        # Define the file path
        file_path_without_dispersion = join(
            data_path_without_dispersion, folder_without_dispersion, "e_crys.txt"
        )
        # Open the file
        with open(file_path_without_dispersion, "r") as file_without_dispersion:
            # Read the value
            value_without_dispersion = float(file_without_dispersion.read())
            # Append the value to the table
            table_without_dispersion.append(
                [folder_without_dispersion, value_without_dispersion]
            )
    df_without_dispersion = pd.DataFrame(
        table_without_dispersion, columns=["folder", "e_crys"]
    )
    df_without_dispersion["e_lattice_without_dispersion"] = (
        df_without_dispersion["e_crys"] - e_gas_without_dispersion
    )
    return (
        data_path_without_dispersion,
        df_without_dispersion,
        e_gas_without_dispersion,
        file_path_without_dispersion,
        file_without_dispersion,
        folder_without_dispersion,
        folders_roman_without_dispersion,
        folders_without_dispersion,
        table_without_dispersion,
        value_without_dispersion,
    )


@app.cell
def __(isdir, join, listdir, pd):
    # Energia con n2p2
    e_gas_n2p2 = -468.46641416393885
    data_path_n2p2 = "../03_crystal_structures_n2p2"
    # List folders
    folders_n2p2 = [
        f for f in listdir(data_path_n2p2) if isdir(join(data_path_n2p2, f))
    ]
    # List folders with roman numerals
    folders_roman_n2p2 = [f for f in folders_n2p2 if f[0].isupper()]
    # Put the value contained in the file e_crys.txt from each directory in a table
    # Define the table
    table_n2p2 = []
    # Loop over the folders
    for folder_n2p2 in folders_roman_n2p2:
        # Define the file path
        file_path_n2p2 = join(data_path_n2p2, folder_n2p2, "e_crys.txt")
        # Open the file
        with open(file_path_n2p2, "r") as file_n2p2:
            # Read the value
            value_n2p2 = float(file_n2p2.read())
            # Append the value to the table
            table_n2p2.append([folder_n2p2, value_n2p2])
    df_n2p2 = pd.DataFrame(table_n2p2, columns=["folder", "e_crys"])
    df_n2p2["e_lattice_n2p2"] = df_n2p2["e_crys"] - e_gas_n2p2
    return (
        data_path_n2p2,
        df_n2p2,
        e_gas_n2p2,
        file_n2p2,
        file_path_n2p2,
        folder_n2p2,
        folders_n2p2,
        folders_roman_n2p2,
        table_n2p2,
        value_n2p2,
    )


@app.cell
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


@app.cell
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


@app.cell
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
    merged["e_lattice_medium"] = merged["e_lattice_medium"] * 1e3
    merged["e_lattice_without_dispersion"] = (
        merged["e_lattice_without_dispersion"] * 1e3
    )
    merged["e_lattice_n2p2"] = merged["e_lattice_n2p2"] * 1e3
    # Calcola il lattice energy in kJ/mol
    merged["e_lattice_dmc_kj"] = merged["e_lattice_dmc"] * conversion_factor
    # Select the row with Ih
    dmc_Ih = merged[merged["folder"] == "Ih"]["e_lattice_dmc_kj"].values[0]
    # Compute the relative lattice energy
    merged["e_lattice_dmc_kj_rel"] = merged["e_lattice_dmc_kj"] - dmc_Ih
    return dmc_Ih, merged


@app.cell
def __(plt, sorted):
    plt.plot(sorted["e_lattice_medium"], label="MACE medium", marker="s", linestyle="-")
    plt.plot(
        sorted["e_lattice_without_dispersion"],
        label="MACE without dispersion",
        marker="^",
        linestyle="-",
    )
    plt.plot(sorted["e_lattice_dmc"], label="DMC", marker="*", linestyle="-")
    plt.plot(sorted["e_lattice_b3lypd4"], label="B3LYP-D4", marker="o", linestyle="-")
    plt.plot(sorted["e_lattice_n2p2"], label="n2p2", marker="x", linestyle="-")
    plt.plot(sorted["e_lattice_lda"], label="LDA", marker="v", linestyle="-")
    plt.xlabel("Crystal structure")
    plt.ylabel("Absolute lattice energy (meV/molecule)")
    plt.xticks(range(len(sorted["folder"])), sorted["folder"])
    plt.legend()
    plt.grid()
    # Make all linewidths smaller
    for line in plt.gca().lines:
        line.set_linewidth(0.5)
    # Make all symbols smaller
    for line in plt.gca().lines:
        line.set_markersize(5)
    plt.title("Absolute lattice energy")
    plt.show()
    return line,


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
def __(conversion_factor, merged):
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

    # Convert e_lattice_medium to kJ/mol
    sorted["e_lattice_medium_kj"] = sorted["e_lattice_medium"] * conversion_factor
    # Select the row with Ih
    medium_Ih = sorted[sorted["folder"] == "Ih"]["e_lattice_medium_kj"].values[0]
    # Compute the relative lattice energy for medium
    sorted["e_lattice_medium_kj_rel"] = sorted["e_lattice_medium_kj"] - medium_Ih

    # Convert e_lattice_without_dispersion to kJ/mol
    sorted["e_lattice_without_dispersion_kj"] = (
        sorted["e_lattice_without_dispersion"] * conversion_factor
    )
    # Select the row with Ih
    without_dispersion_Ih = sorted[sorted["folder"] == "Ih"][
        "e_lattice_without_dispersion_kj"
    ].values[0]
    # Compute the relative lattice energy for without dispersion
    sorted["e_lattice_without_dispersion_kj_rel"] = (
        sorted["e_lattice_without_dispersion_kj"] - without_dispersion_Ih
    )

    # Convert n2p2 to kJ/mol
    sorted["e_lattice_n2p2_kj"] = sorted["e_lattice_n2p2"] * conversion_factor
    # Select the row with Ih
    n2p2_Ih = sorted[sorted["folder"] == "Ih"]["e_lattice_n2p2_kj"].values[0]
    # Compute the relative lattice energy for n2p2
    sorted["e_lattice_n2p2_kj_rel"] = sorted["e_lattice_n2p2_kj"] - n2p2_Ih
    return medium_Ih, n2p2_Ih, sorted, sorting_order, without_dispersion_Ih


@app.cell
def __(plt, sorted):
    # Plot the relative lattice energy
    plt.plot(
        sorted["e_lattice_dmc_kj_rel"], label="DMC", marker="o", linestyle="-"
    )
    plt.plot(
        sorted["e_lattice_medium_kj_rel"],
        label="MACE medium",
        marker="s",
        linestyle="-",
    )
    plt.plot(
        sorted["e_lattice_without_dispersion_kj_rel"],
        label="MACE without dispersion",
        marker="^",
        linestyle="-",
    )
    plt.plot(
        sorted["e_lattice_n2p2_kj_rel"],
        label="n2p2",
        marker="x",
        linestyle="-",
    )
    plt.xlabel("Crystal structure")
    plt.ylabel("Relative lattice energy (kJ/mol)")
    plt.xticks(range(len(sorted["folder"])), sorted["folder"])
    plt.legend()
    plt.grid()
    plt.title("Relative lattice energy")
    plt.show()
    return


@app.cell
def __(plt, sorted):
    plt.plot(sorted["e_lattice_medium_kj"], label="MACE medium", marker="s", linestyle="-")
    plt.plot(
        sorted["e_lattice_without_dispersion_kj"],
        label="MACE without dispersion",
        marker="^",
        linestyle="-",
    )
    plt.plot(sorted["e_lattice_dmc_kj"], label="DMC", marker="*", linestyle="-")
    # plt.plot(sorted["e_lattice_b3lypd4_kj"], label="B3LYP-D4", marker="o", linestyle="-")
    plt.plot(sorted["e_lattice_n2p2_kj"], label="n2p2", marker="x", linestyle="-")
    # plt.plot(sorted["e_lattice_lda_kj"], label="LDA", marker="v", linestyle="-")
    plt.xlabel("Crystal structure")
    plt.ylabel("Absolute lattice energy (meV/molecule)")
    plt.xticks(range(len(sorted["folder"])), sorted["folder"])
    plt.legend()
    plt.grid()
    plt.title("Absolute lattice energy")
    plt.show()
    return


if __name__ == "__main__":
    app.run()
