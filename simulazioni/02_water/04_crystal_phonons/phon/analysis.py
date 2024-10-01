import marimo

__generated_with = "0.8.22"
app = marimo.App(width="medium", layout_file="layouts/analysis.grid.json")


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
def __(mo, pd, plt):
    def plot_bandpath(filepath):
        ticks = [0.00000, 0.06911, 0.18009, 0.24920, 0.33082, 0.39993, 0.50206]
        labels_gupta = ["G", "A", "K", "H", "M", "L", "G"]
        freq = pd.read_csv(
            filepath, sep="[ ]+", header=None, engine="python", index_col=0
        )
        plt.plot(freq, color="green")
        plt.ylim(-0.5, 5)
        ax = plt.gca()
        ax.set_xticks(ticks, labels=labels_gupta)
        plt.grid(axis="x")
        plt.ylabel("Frequency (THz)")
        return mo.mpl.interactive(plt.gcf())
    return (plot_bandpath,)


@app.cell
def __(file, plot_bandpath):
    plot_bandpath(file.path())
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
def __(file_ref, plot_bandpath):
    plot_bandpath(file_ref.path())
    return


@app.cell
def __(mo, pd, plt):
    from scipy.ndimage import gaussian_filter1d

    def plot_dos(filepath, sigma=1.0):
        dos = pd.read_csv(
            filepath, sep="[ ]+", header=None, engine="python", index_col=0
        )
        dos_smeared = gaussian_filter1d(dos.iloc[:, 0], sigma=sigma)
        dos_smeared_df = pd.DataFrame({"Energy": dos.index, "DOS": dos_smeared})
        plt.plot(dos_smeared_df["Energy"], dos_smeared_df["DOS"], color="blue")
        plt.ylabel("DOS")
        plt.xlim(0, 11)
        plt.xlabel("Frequency (THz)")
        plt.grid(axis="x")
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


if __name__ == "__main__":
    app.run()
