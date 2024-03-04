import marimo

__generated_with = "0.2.13"
app = marimo.App()


@app.cell
def __():
    import marimo as mo
    import pandas as pd
    import matplotlib.pyplot as plt
    import glob
    return glob, mo, pd, plt


@app.cell
def __(mo):
    mo.md("# Ottimizzazione della geometria del dimero d'acqua e studio della stabilità")
    return


@app.cell
def __(mo):
    models = ["small", "medium", "large"]
    model = mo.ui.dropdown(options=models, value="small", label="Modello")
    model
    return model, models


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
def __(glob):
    def get_filenames(model):
        root = "MACE-MP-0"
        filenames = glob.glob(f"{root}/{model}/*_summary.txt")
        return sorted(
            filenames, key=lambda x: float(x.split("delta=")[1].split("_")[0])
        )
    return get_filenames,


@app.cell
def __(immaginari_a_negativi, pd):
    def get_frequencies(filename):
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
    return get_frequencies,


@app.cell
def __(get_filenames, model):
    filenames = get_filenames(model.value)
    return filenames,


@app.cell
def __(get_frequencies):
    def build_frequencies_dataframe(filename):
        df = get_frequencies(filename)
        df["delta"] = float(filename.split("delta=")[1].split("_")[0])
        return df
    return build_frequencies_dataframe,


@app.cell
def __(build_frequencies_dataframe, filenames, pd):
    df = pd.DataFrame()
    for _f in filenames:
        _df = build_frequencies_dataframe(_f)
        df = pd.concat([df, _df])
    return df,


@app.cell
def __(df, model, plt):
    groups = df.groupby("#")
    for _name, _group in groups:
        plt.plot(_group["delta"], _group["cm^-1"], marker="x", label=_name)
    plt.xscale("log")
    plt.yscale("symlog", linthresh=1e-1)
    plt.legend(ncol=2)
    plt.ylim(-1e4, 1e4)
    plt.title(f"MACE-MP-0 {model.value} D")
    plt.xlabel("Displacement (Å)")
    plt.ylabel("Frequency (cm^-1)")
    return groups,


if __name__ == "__main__":
    app.run()
