import marimo

__generated_with = "0.2.13"
app = marimo.App(layout_file="layouts/01.9_analysis.grid.json")


@app.cell
def __():
    import marimo as mo
    import pandas as pd
    import matplotlib.pyplot as plt
    import glob
    return glob, mo, pd, plt


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        """
        # Analisi dei risultati

        Ho avviato run con e senza dispersione, per i tre modelli small, medium e large, e per valori diversi della fmax di convergenza. Dopodiché ho studiato le vibrazioni con diversi displacement per ciascun caso.

        Iniziamo analizzando un caso campione rispetto al displacement per il calcolo delle vibrazioni.
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
        if dispersion:
            root = "01.8_molecule_converge_fmax_parallel_dispersion"
        else:
            root = "01.7_molecule_converge_fmax_parallel"
        filenames = glob.glob(f"{root}/" f"{model}/{fmax}/" "*_summary.txt")
        filenames = sorted(
            filenames, key=lambda x: float(x.split("delta=")[1].split("_")[0])
        )
        return filenames
    return get_summaries,


@app.cell
def __(mo):
    models = ["small", "medium", "large"]
    model = mo.ui.dropdown(options=models, value="small", label="Modello")

    fmaxs = ["1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"]
    fmax = mo.ui.dropdown(options=fmaxs, value="1e-4", label="fmax")

    dispersion = mo.ui.checkbox(value=False, label="Dispersione")

    mo.hstack([model, fmax, dispersion])
    return dispersion, fmax, fmaxs, model, models


@app.cell(hide_code=True)
def __(dispersion, fmax, get_summaries, model, pd, read_frequencies):
    filenames = get_summaries(model.value, fmax.value, dispersion.value)

    df = pd.DataFrame()
    for _f in filenames:
        _df = read_frequencies(_f)
        _df["delta"] = float(_f.split("delta=")[1].split("_")[0])
        df = pd.concat([df, _df])
    return df, filenames


@app.cell
def __(df, plt):
    groups = df.groupby("#")
    for _name, _group in groups:
        plt.plot(_group["delta"], _group["cm^-1"], marker="x", label=_name)
    plt.xscale("log")
    plt.yscale("symlog", linthresh=1e-1)
    plt.legend(ncol=2)
    return groups,


if __name__ == "__main__":
    app.run()
