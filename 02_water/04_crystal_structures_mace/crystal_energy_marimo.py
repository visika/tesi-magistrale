import marimo

__generated_with = "0.2.5"
app = marimo.App()


@app.cell
def __():
    # Import listdir
    from os import listdir
    from os.path import isfile, join, isdir

    import matplotlib.pyplot as plt
    return isdir, isfile, join, listdir, plt


@app.cell
def __():
    e_gas = -14.160290476990022
    e_gas_large = -14.15942157689177
    e_gas_n2p2 = -468.46641416393885
    return e_gas, e_gas_large, e_gas_n2p2


@app.cell
def __(isdir, join, listdir):
    # Define data_path
    data_path = 'medium'
    # List folders
    folders = [f for f in listdir(data_path) if isdir(join(data_path, f))]
    # List folders with roman numerals
    folders_roman = [f for f in folders if f[0].isupper()]
    return data_path, folders, folders_roman


@app.cell
def __(data_path, folders_roman, join):
    # Put the value contained in the file e_crys.txt from each directory in a table
    # Define the table
    table = []
    # Loop over the folders
    for folder in folders_roman:
        # Define the file path
        file_path = join(data_path, folder, 'e_crys.txt')
        # Open the file
        with open(file_path, 'r') as file:
            # Read the value
            value = float(file.read())
            # Append the value to the table
            table.append([folder, value])
    return file, file_path, folder, table, value


@app.cell
def __(table):
    # Render the table as a pandas dataframe
    import pandas as pd
    df = pd.DataFrame(table, columns=['folder', 'e_crys'])
    return df, pd


@app.cell
def __(e_gas_large, isdir, join, listdir, pd):
    # Define data_path
    data_path_large = "large"
    # List folders
    folders_large = [f for f in listdir(data_path_large) if isdir(join(data_path_large, f))]
    # List folders with roman numerals
    folders_roman_large = [f for f in folders_large if f[0].isupper()]
    # Put the value contained in the file e_crys.txt from each directory in a table
    # Define the table
    table_large = []
    # Loop over the folders
    for folder_large in folders_roman_large:
        # Define the file path
        file_path_large = join(data_path_large, folder_large, 'e_crys.txt')
        # Open the file
        with open(file_path_large, 'r') as file_large:
            # Read the value
            value_large = float(file_large.read())
            # Append the value to the table
            table_large.append([folder_large, value_large])
    df_large = pd.DataFrame(table_large, columns=['folder', 'e_crys'])
    df_large["e_lattice_large"] = df_large["e_crys"] - e_gas_large
    return (
        data_path_large,
        df_large,
        file_large,
        file_path_large,
        folder_large,
        folders_large,
        folders_roman_large,
        table_large,
        value_large,
    )


@app.cell
def __(e_gas_n2p2, isdir, join, listdir, pd):
    data_path_n2p2 = "../03_crystal_structures_n2p2"
    # List folders
    folders_n2p2 = [f for f in listdir(data_path_n2p2) if isdir(join(data_path_n2p2, f))]
    # List folders with roman numerals
    folders_roman_n2p2 = [f for f in folders_n2p2 if f[0].isupper()]
    # Put the value contained in the file e_crys.txt from each directory in a table
    # Define the table
    table_n2p2 = []
    # Loop over the folders
    for folder_n2p2 in folders_roman_n2p2:
        # Define the file path
        file_path_n2p2 = join(data_path_n2p2, folder_n2p2, 'e_crys.txt')
        # Open the file
        with open(file_path_n2p2, 'r') as file_n2p2:
            # Read the value
            value_n2p2 = float(file_n2p2.read())
            # Append the value to the table
            table_n2p2.append([folder_n2p2, value_n2p2])
    df_n2p2 = pd.DataFrame(table_n2p2, columns=['folder', 'e_crys'])
    df_n2p2["e_lattice_n2p2"] = df_n2p2["e_crys"] - e_gas_n2p2
    return (
        data_path_n2p2,
        df_n2p2,
        file_n2p2,
        file_path_n2p2,
        folder_n2p2,
        folders_n2p2,
        folders_roman_n2p2,
        table_n2p2,
        value_n2p2,
    )


@app.cell
def __(df):
    # Calculate the difference of the values in the column e_crys with respect to the value associated to the folder Ih
    # Define the reference value
    ref_value = df[df['folder'] == 'Ih']['e_crys'].values[0]
    # Calculate the difference
    df['diff'] = df['e_crys'] - ref_value
    return ref_value,


@app.cell
def __(df, e_gas):
    # Calculate the lattice energy as the difference between e_crys and e_gas
    df['e_lattice'] = df['e_crys'] - e_gas
    return


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
    df_b3lypd4 = pd.DataFrame(
        b3lypd4_table, columns=["folder", "e_lattice_b3lypd4"]
    )
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
    # Expected: 1 meV/particle = 0.0965 kJ/mol
    dmc_kj = -59.45
    dmc_ev = -616
    conversion_factor = dmc_kj / dmc_ev
    return conversion_factor, dmc_ev, dmc_kj


@app.cell
def __(
    conversion_factor,
    df,
    df_b3lypd4,
    df_dmc,
    df_large,
    df_lda,
    df_n2p2,
    pd,
):
    merged = pd.merge(df, df_dmc, on="folder", suffixes=("_medium", "_dmc"))
    merged = pd.merge(merged, df_b3lypd4, on="folder")
    merged = pd.merge(
        merged, df_large, on="folder", suffixes=("_merged", "_large")
    )
    merged = pd.merge(merged, df_n2p2, on="folder", suffixes=("_merged", "_n2p2"))
    merged = pd.merge(merged, df_lda, on="folder", suffixes=("_merged", "_lda"))
    merged["e_lattice_medium"] = merged["e_lattice_medium"] * 1e3
    merged["e_lattice_large"] = merged["e_lattice_large"] * 1e3
    merged["e_lattice_n2p2"] = merged["e_lattice_n2p2"] * 1e3

    # Calcola il lattice energy in kJ/mol
    merged["e_lattice_dmc_kj"] = merged["e_lattice_dmc"] * conversion_factor

    # Select the row with Ih
    dmc_Ih = merged[merged["folder"] == "Ih"]["e_lattice_dmc_kj"].values[0]
    # Compute the relative lattice energy
    merged["e_lattice_dmc_kj_rel"] = merged["e_lattice_dmc_kj"] - dmc_Ih

    # Calcola la lattice energy di B3LYP-D4 in kJ/mol
    merged["e_lattice_b3lypd4_kj"] = merged["e_lattice_b3lypd4"] * conversion_factor
    return dmc_Ih, merged


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

    # Convert e_lattice_large to kJ/mol
    sorted["e_lattice_large_kj"] = sorted["e_lattice_large"] * conversion_factor
    # Select the row with Ih
    large_Ih = sorted[sorted["folder"] == "Ih"]["e_lattice_large_kj"].values[0]
    # Compute the relative lattice energy for large
    sorted["e_lattice_large_kj_rel"] = sorted["e_lattice_large_kj"] - large_Ih

    # Convert n2p2 to kJ/mol
    sorted["e_lattice_n2p2_kj"] = sorted["e_lattice_n2p2"] * conversion_factor
    # Select the row with Ih
    n2p2_Ih = sorted[sorted["folder"] == "Ih"]["e_lattice_n2p2_kj"].values[0]
    # Compute the relative lattice energy for n2p2
    sorted["e_lattice_n2p2_kj_rel"] = sorted["e_lattice_n2p2_kj"] - n2p2_Ih

    # Compute e_lattice_b3lypd4_kj_rel
    b3lypd4_Ih = sorted[sorted["folder"] == "Ih"]["e_lattice_b3lypd4_kj"].values[0]
    sorted["e_lattice_b3lypd4_kj_rel"] = sorted["e_lattice_b3lypd4_kj"] - b3lypd4_Ih

    # Convert LDAs to kJ/mol
    sorted["e_lattice_lda_kj"] = sorted["e_lattice_lda"] * conversion_factor
    # Select the row with Ih
    lda_Ih = sorted[sorted["folder"] == "Ih"]["e_lattice_lda_kj"].values[0]
    # Compute the relative lattice energy for LDA
    sorted["e_lattice_lda_kj_rel"] = sorted["e_lattice_lda_kj"] - lda_Ih
    return (
        b3lypd4_Ih,
        large_Ih,
        lda_Ih,
        medium_Ih,
        n2p2_Ih,
        sorted,
        sorting_order,
    )


@app.cell(hide_code=True)
def __(plt, sorted):
    plt.plot(
        sorted["e_lattice_dmc"],
        label="DMC",
        marker="*",
        linestyle="-",
        markersize=10,
        linewidth=0.7,
        color="#FFDE00",
    )
    markersize = 5
    linewidth = 0.5
    plt.plot(
        sorted["e_lattice_b3lypd4"],
        label="B3LYP-D4",
        marker="o",
        linestyle="-",
        markersize=markersize,
        linewidth=linewidth,
    )
    plt.plot(
        sorted["e_lattice_medium"],
        label="MACE medium",
        marker="s",
        linestyle="-",
        markersize=markersize,
        linewidth=linewidth,
    )
    plt.plot(
        sorted["e_lattice_large"],
        label="MACE large",
        marker="^",
        linestyle="-",
        markersize=markersize,
        linewidth=linewidth,
    )
    plt.plot(
        sorted["e_lattice_n2p2"],
        label="n2p2",
        marker="x",
        linestyle="-",
        markersize=markersize,
        linewidth=linewidth,
    )
    plt.plot(
        sorted["e_lattice_lda"],
        label="LDA",
        marker="v",
        linestyle="-",
        markersize=markersize,
        linewidth=linewidth,
    )
    plt.xlabel("Crystal structure")
    plt.ylabel("Absolute lattice energy (meV/molecule)")
    plt.xticks(range(len(sorted["folder"])), sorted["folder"])
    plt.legend()
    plt.grid()
    plt.title("Absolute lattice energy")
    return linewidth, markersize


@app.cell(hide_code=True)
def __(sorted):
    # Compute the mean absolute error in kJ/mol
    mae_medium_kj = (
        (sorted["e_lattice_medium_kj"] - sorted["e_lattice_dmc_kj"]).abs().mean()
    )
    mae_large_kj = (
        (sorted["e_lattice_large_kj"] - sorted["e_lattice_dmc_kj"]).abs().mean()
    )
    mae_n2p2_kj = (
        (sorted["e_lattice_n2p2_kj"] - sorted["e_lattice_dmc_kj"]).abs().mean()
    )
    mae_b3lypd4_kj = (
        (sorted["e_lattice_b3lypd4_kj"] - sorted["e_lattice_dmc_kj"]).abs().mean()
    )

    import marimo as mo

    mo.md(
        f"""
        ## MAE for the absolute lattice energy
        MACE medium: {mae_medium_kj.round(2)} kJ/mol  
        MACE large: {mae_large_kj.round(2)} kJ/mol  
        n2p2: {mae_n2p2_kj.round(2)} kJ/mol  
        B3LYP-D4: {mae_b3lypd4_kj.round(2)} kJ/mol
    """
    )
    return mae_b3lypd4_kj, mae_large_kj, mae_medium_kj, mae_n2p2_kj, mo


@app.cell(hide_code=True)
def __(merged):
    # Compute the mean absolute error in eV/molecule
    mae_medium = (
        (merged["e_lattice_medium"] - merged["e_lattice_dmc"]).abs().mean()
    )
    mae_large = (merged["e_lattice_large"] - merged["e_lattice_dmc"]).abs().mean()
    mae_b3lypd4 = (
        (merged["e_lattice_b3lypd4"] - merged["e_lattice_dmc"]).abs().mean()
    )
    mae_n2p2 = (merged["e_lattice_n2p2"] - merged["e_lattice_dmc"]).abs().mean()
    # mae_medium, mae_large, mae_b3lypd4, mae_n2p2
    return mae_b3lypd4, mae_large, mae_medium, mae_n2p2


@app.cell(hide_code=True)
def __(mae_b3lypd4, mae_large, mae_medium, mae_n2p2):
    def get_variable_name(variable):
        variable_name = [
            name for name, value in globals().items() if value is variable
        ][0]
        return variable_name


    for e in [mae_medium, mae_large, mae_b3lypd4, mae_n2p2]:
        print(f"{get_variable_name(e)}: {e}")
    return e, get_variable_name


@app.cell(hide_code=True)
def __(plt, sorted):
    # Plot the relative lattice energy
    plt.plot(
        sorted["e_lattice_dmc_kj_rel"],
        label="DMC",
        marker="*",
        linestyle="-",
        markersize=10,
        linewidth=0.7,
        color="#FFDE00",
    )
    plt.plot(
        sorted["e_lattice_b3lypd4_kj_rel"],
        label="B3LYP-D4",
        marker="o",
        linestyle="-",
        markersize=5,
        linewidth=0.5,
    )
    plt.plot(
        sorted["e_lattice_medium_kj_rel"],
        label="MACE medium",
        marker="s",
        linestyle="-",
        markersize=5,
        linewidth=0.5,
    )
    plt.plot(
        sorted["e_lattice_large_kj_rel"],
        label="MACE large",
        marker="^",
        linestyle="-",
        markersize=5,
        linewidth=0.5,
    )
    plt.plot(
        sorted["e_lattice_n2p2_kj_rel"],
        label="n2p2",
        marker="x",
        linestyle="-",
        markersize=5,
        linewidth=0.5,
    )
    plt.plot(
        sorted["e_lattice_lda_kj_rel"],
        label="LDA",
        marker="v",
        linestyle="-",
        markersize=5,
        linewidth=0.5,
    )
    plt.xlabel("Crystal structure")
    plt.ylabel("Relative lattice energy (kJ/mol)")
    plt.xticks(range(len(sorted["folder"])), sorted["folder"])
    plt.legend()
    plt.grid()
    plt.title("Relative lattice energy")
    return


@app.cell(hide_code=True)
def __(mo, sorted):
    # Compute the MAE for the relative lattice energy
    mae_medium_kj_rel = (
        (sorted["e_lattice_medium_kj_rel"] - sorted["e_lattice_dmc_kj_rel"])
        .abs()
        .mean()
    )
    mae_large_kj_rel = (
        (sorted["e_lattice_large_kj_rel"] - sorted["e_lattice_dmc_kj_rel"])
        .abs()
        .mean()
    )
    mae_n2p2_kj_rel = (
        (sorted["e_lattice_n2p2_kj_rel"] - sorted["e_lattice_dmc_kj_rel"])
        .abs()
        .mean()
    )
    mae_b3lypd4_kj_rel = (
        (sorted["e_lattice_b3lypd4_kj_rel"] - sorted["e_lattice_dmc_kj_rel"])
        .abs()
        .mean()
    )
    mo.md(
        f"""
        ## MAE for the relative lattice energy
        MACE medium: {mae_medium_kj_rel.round(2)} kJ/mol  
        MACE large: {mae_large_kj_rel.round(2)} kJ/mol  
        n2p2: {mae_n2p2_kj_rel.round(2)} kJ/mol  
        B3LYP-D4: {mae_b3lypd4_kj_rel.round(2)} kJ/mol
    """
    )
    return (
        mae_b3lypd4_kj_rel,
        mae_large_kj_rel,
        mae_medium_kj_rel,
        mae_n2p2_kj_rel,
    )


if __name__ == "__main__":
    app.run()
