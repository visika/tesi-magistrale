import marimo

__generated_with = "0.3.3"
app = marimo.App(layout_file="layouts/01.9_analysis.grid.json")


@app.cell
def __():
    import marimo as mo
    import pandas as pd
    import matplotlib.pyplot as plt
    import glob
    return glob, mo, pd, plt


@app.cell
def __(mo):
    mo.md(
        """
        # Analisi dei risultati sulla molecola d'acqua

        Ho avviato run di ottimizzazione della molecola di \(\mathrm{H_2O}\) con e senza dispersione, per i tre modelli small, medium e large di MACE-MP-0, e per valori diversi della `fmax` di convergenza dell'ottimizzatore. Dopodiché ho studiato le vibrazioni con diversi displacement per ciascun caso.
        """
    )
    return


@app.cell
def __():
    def immaginari_a_negativi(stringa):
        """
        Se c'è l'unità immaginaria, restituisci il numero reale negativo.
        Se non c'è l'unità immaginaria, restituisci il numero reale positivo.
        """
        if "i" in str(stringa):
            return -float(stringa.split("i")[0])
        else:
            return float(stringa)
    return immaginari_a_negativi,


@app.cell
def __(immaginari_a_negativi, pd):
    def read_frequencies(filename):
        with open(
            filename,
            "r",
        ) as _f:
            df = pd.read_csv(
                _f.name,
                skiprows=[0, 2],
                skipfooter=2,
                engine="python",
                delimiter="\ +",
            )
        df["cm^-1"] = df["cm^-1"].apply(immaginari_a_negativi)
        return df
    return read_frequencies,


@app.cell
def __(glob):
    def get_summaries(model, fmax, dispersion):
        if model in ["small", "medium", "large"]:
            if dispersion:
                root = "01.8_molecule_converge_fmax_parallel_dispersion"
            else:
                root = "01.7_molecule_converge_fmax_parallel"
            filenames = glob.glob(f"{root}/" f"{model}/{fmax}/" "*_summary.txt")
            filenames = sorted(
                filenames, key=lambda x: float(x.split("delta=")[1].split("_")[0])
            )
        elif model == "MACE-ICE13":
            root = model
            filenames = glob.glob(
                f"{root}/" f"dispersion={dispersion}/" f"{fmax}/*_summary.txt"
            )
            filenames = sorted(
                filenames, key=lambda x: float(x.split("delta=")[1].split("_")[0])
            )
        elif model == "MACE-ICE13-1":
            root = model
            filenames = glob.glob(f"{root}/" f"*_summary.txt")
            filenames = sorted(
                filenames, key=lambda x: float(x.split("delta=")[1].split("_")[0])
            )
        return filenames
    return get_summaries,


@app.cell
def __(mo):
    models = ["small", "medium", "large", "MACE-ICE13", "MACE-ICE13-1"]
    model = mo.ui.dropdown(options=models, value="small", label="Modello")

    fmaxs = ["1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"]
    fmax = mo.ui.dropdown(options=fmaxs, value="1e-4", label="fmax")

    dispersion = mo.ui.checkbox(value=False, label="Dispersione")

    mo.hstack([model, fmax, dispersion])
    return dispersion, fmax, fmaxs, model, models


@app.cell
def __(dispersion, fmax, get_summaries, model, pd, read_frequencies):
    filenames = get_summaries(model.value, fmax.value, dispersion.value)

    df = pd.DataFrame()
    for _f in filenames:
        _df = read_frequencies(_f)
        _df["delta"] = float(_f.split("delta=")[1].split("_")[0])
        df = pd.concat([df, _df])
    return df, filenames


@app.cell
def __(df, dispersion, fmax, model, plt):
    groups = df.groupby("#")
    for _name, _group in groups:
        plt.plot(_group["delta"], _group["cm^-1"], marker="x", label=_name)
    plt.xscale("log")
    plt.yscale("symlog", linthresh=1e-1)
    plt.ylim(-1e4, 1e4)
    plt.legend(ncol=2)

    if model.value == "MACE-ICE13-1":
        title_fmax = ", fmax=1e-8"
    else:
        title_fmax = f", fmax={fmax.value}"

    plt.title(
        f"{'MACE-MP-0' if model.value in ['small', 'medium', 'large'] else ''} {model.value}"
        f"{' D' if dispersion.value else ''}"
        f"{title_fmax}"
    )
    plt.xlabel("Displacement (Å)")
    plt.ylabel("Frequency (cm^-1)")

    # plt.savefig("MACE-xxx-y/frequencies.png")
    return groups, title_fmax


if __name__ == "__main__":
    app.run()
