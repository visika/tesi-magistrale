import marimo

__generated_with = "0.1.81"
app = marimo.App()


@app.cell
def __():
    import pandas as pd
    return pd,


@app.cell
def __():
    import glob

    # List filenames ending in .txt
    filenames = glob.glob("*.txt")
    filenames
    # Sort the filenames based on the number
    filenames_sorted = sorted(
        filenames, key=lambda x: float(x.split("=")[1].split("_")[0])
    )
    filenames_sorted
    return filenames, filenames_sorted, glob


@app.cell
def __(filenames_sorted, immagina, pd):
    df = pd.DataFrame()
    for filename in filenames_sorted:
        # Read file into a DataFrame: df
        df_tmp = pd.read_csv(
            filename,
            sep="\ +",
            header=0,
            skiprows=[0, 2],
            skipfooter=2,
            engine="python",
        )
        df_tmp["delta"] = float(filename.split("=")[1].split("_")[0])
        df = pd.concat([df, df_tmp])

    frequenze = [immagina(stringa) for stringa in df["cm^-1"].values]
    # Replace column cm^-1 with frequenze
    df["freq"] = frequenze
    return df, df_tmp, filename, frequenze


@app.cell
def __(df):
    # Plot frequencies vs delta grouped by number
    import matplotlib.pyplot as plt

    # Define size
    plt.figure(figsize=(7, 5))
    groups = df.groupby("#")
    for name, group in groups:
        plt.plot(
            group["delta"], group["freq"], marker="x", linestyle="-", label=name
        )
    # Logscale on x
    plt.xscale("log")
    # Symmetrical logscale on y
    plt.yscale("symlog", linthresh=0.1)
    # Legend with two columns
    plt.legend(ncol=2)
    # Grid
    plt.grid()
    plt.ylim(bottom=-1e4, top=1e4)
    # Title
    plt.title(
        "MACE-MP-0 small, water dimer: frequencies vs displacement\n"
        "keeping the positions in memory"
    )
    # Labels
    # delta in angstrom with latex
    plt.xlabel("$\delta$ ($\AA$)")
    plt.ylabel("Frequenze ($\mathrm{cm}^{-1}$)")
    # Save figure
    # plt.savefig("frequenze_vs_delta.png", dpi=300)
    # Show figure
    plt.show()
    return group, groups, name, plt


@app.cell
def __():
    def immagina(stringa):
        """
        Se c'è l'unità immaginaria, restituisci il numero reale negativo.
        Se non c'è l'unità immaginaria, restituisci il numero reale positivo.
        """
        if "i" in stringa:
            return -float(stringa.split("i")[0])
        else:
            return float(stringa)
    return immagina,


if __name__ == "__main__":
    app.run()
