import marimo

__generated_with = "0.8.22"
app = marimo.App(width="full", layout_file="layouts/analysis.grid.json")


@app.cell
def __():
    # Matplotlib per i grafici
    import matplotlib.pyplot as plt
    # Pandas per organizzare i dati
    import pandas as pd
    # Marimo per l'interfaccia utente
    import marimo as mo
    return mo, pd, plt


@app.cell
def __(mo):
    mo.md(r"""## Struttura a bande calcolata con MACE""")
    return


@app.cell
def __(mo):
    # Default: 16.MACE_geometrie_Flaviano/FREQ1
    file = mo.ui.file_browser()
    file
    return (file,)


@app.cell
def __(pd, plt):
    def plot_bandpath(filepath, title=""):
        ticks = [0.00000, 0.06911, 0.18009, 0.24920, 0.33082, 0.39993, 0.50206]
        labels_gupta = ["G", "A", "K", "H", "M", "L", "G"]
        freq = pd.read_csv(
            filepath, sep="[ ]+", header=None, engine="python", index_col=0
        )
        plt.plot(freq, color="green")
        plt.ylim(-0.5, 11)
        ax = plt.gca()
        ax.set_xticks(ticks, labels=labels_gupta)
        plt.grid(axis="x")
        plt.ylabel("Frequenza (THz)")
        plt.suptitle(title)
        return plt.gcf()
    return (plot_bandpath,)


@app.cell
def __(mo, plot_bandpath):
    # mo.mpl.interactive(plot_bandpath(file.path()))
    mo.mpl.interactive(plot_bandpath("16.MACE_geometrie_Flaviano/FREQ1"))
    return


@app.cell
def __(mo):
    mo.md(r"""## Struttura a bande calcolata con revPBE-D3""")
    return


@app.cell
def __(mo):
    # Default: 15.revPBED3/FREQ1
    file_ref = mo.ui.file_browser()
    file_ref
    return (file_ref,)


@app.cell
def __(mo, plot_bandpath):
    # plot_bandpath(file_ref.path())
    mo.mpl.interactive(plot_bandpath("15.revPBED3/FREQ1"))
    return


@app.cell
def __(mo, pd, plt):
    from scipy.ndimage import gaussian_filter1d


    def plot_dos(filepath, sigma=1.0):
        # Set figure dimensions
        plt.figure(figsize=(6, 1))
        
        dos = pd.read_csv(
            filepath, sep="[ ]+", header=None, engine="python", index_col=0
        )
        dos_smeared = gaussian_filter1d(dos.iloc[:, 0], sigma=sigma)
        dos_smeared_df = pd.DataFrame({"Energy": dos.index, "DOS": dos_smeared})
        # plt.plot(dos_smeared_df["Energy"], dos_smeared_df["DOS"], color="blue")
        plt.fill_between(
            dos_smeared_df["Energy"],
            dos_smeared_df["DOS"],
            y2=0,
            color="grey",
            edgecolor="k",
            lw=1,
        )
        # plt.ylabel("DOS")
        plt.xlim(-0.5, 11)
        plt.ylim(bottom=0)
        # plt.xlabel("Frequency (THz)")
        plt.xticks([])
        plt.yticks([])
        # plt.grid(axis="x")
        return mo.mpl.interactive(plt.gcf())
    return gaussian_filter1d, plot_dos


@app.cell
def __(mo):
    mo.md(r"""## DOS calcolata con MACE""")
    return


@app.cell
def __(plot_dos):
    plot_dos("16.MACE_geometrie_Flaviano/DOS")
    return


@app.cell
def __(mo):
    mo.md(r"""## DOS calcolata con revPBE-D3""")
    return


@app.cell
def __(plot_dos):
    plot_dos("15.revPBED3/DOS")
    return


@app.cell
def __(mo):
    mo.md(r"""## Grafici dei diagrammi a bande con la DOS a lato""")
    return


@app.cell
def __(ax2, gaussian_filter1d, pd, plt, transforms):
    def plot_band_dos(directory, max_freq):
        # fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
        fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))

        ticks = [0.00000, 0.06911, 0.18009, 0.24920, 0.33082, 0.39993, 0.50206]
        labels_gupta = ["G", "A", "K", "H", "M", "L", "G"]
        freq = pd.read_csv(
            f"{directory}/FREQ1",
            sep="[ ]+",
            header=None,
            engine="python",
            index_col=0,
        )
        ax1.plot(freq, color="green")
        ax1.set_ylim(-0.5, 5)
        ax1.set_xticks(ticks, labels=labels_gupta)
        ax1.grid(axis="x")
        ax1.set_ylabel("Frequency (THz)")

        dos = pd.read_csv(
            f"{directory}/DOS",
            sep="[ ]+",
            header=None,
            engine="python",
            index_col=0,
        )
        dos_smeared = gaussian_filter1d(dos.iloc[:, 0], sigma=1.0)
        dos_smeared_df = pd.DataFrame({"Energy": dos.index, "DOS": dos_smeared})
        # Sort the values of the dataframe
        dos_smeared_df = dos_smeared_df.sort_values(by="Energy")

        # ax2.plot(dos_smeared_df["Energy"], dos_smeared_df["DOS"], color="blue")
        base = ax2.transData
        rot = transforms.Affine2D().rotate_deg(90)

        ax2.fill_between(
            dos_smeared_df["Energy"],
            dos_smeared_df["DOS"],
            y2=0,
            color="grey",
            edgecolor="k",
            lw=1,
            transform=rot + base,
        )
        ax2.set_ylabel("Energy")
        # ax2.set_xlim(left=0)
        # ax2.set_ylim(-0.5, 5)
        # ax2.set_yticks([])
        # ax2.set_xticks([])
        ax2.set_xlabel("Frequency (THz)")
        # ax2.grid(axis="x")


        return plt.gcf()
    return (plot_band_dos,)


if __name__ == "__main__":
    app.run()
