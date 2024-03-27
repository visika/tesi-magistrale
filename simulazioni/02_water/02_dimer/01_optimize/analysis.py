import marimo

__generated_with = "0.3.3"
app = marimo.App(layout_file="layouts/analysis.grid.json")


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
    models = ["small", "medium", "large", "MACE-ICE13", "MACE-ICE13-1", "n2p2"]
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
        if model in ["small", "medium", "large"]:
            root = "MACE-MP-0"
            filenames = glob.glob(f"{root}/{model}/*_summary.txt")
        elif model in ["MACE-ICE13", "MACE-ICE13-1"]:
            root = model
            filenames = glob.glob(f"{root}/*_summary.txt")
        elif model == "n2p2":
            root = "n2p2"
            filenames = glob.glob(f"{root}/*_summary.txt")
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
    if model.value in ["small", "medium", "large"]:
        plt.title(f"MACE-MP-0 {model.value} D")
    elif model.value in ["MACE-ICE13", "MACE-ICE13-1"]:
        plt.title(model.value)
    elif model.value == "n2p2":
        plt.title(f"n2p2")
    plt.xlabel("Displacement (Å)")
    plt.ylabel("Frequency (cm^-1)")
    return groups,


@app.cell(hide_code=True)
def __(mo):
    mo.md(
        """
        ## Studio della geometria

        Ottieni gli indici selezionando la voce di menu:
        `View > Show Labels > Atom Index`

        - 0 ossigeno accettore
        - 1 idrogeno accettore
        - 2 idrogeno accettore
        - 3 ossigeno donore
        - 4 idrogeno donore esterno (f)
        - 5 idrogeno donore che punta (d)
        """
    )
    return


@app.cell
def __():
    from ase.io import read
    from ase.visualize import view
    return read, view


@app.cell
def __(model, read):
    def atoms_reader(model: str):
        if model in ["small", "medium", "large"]:
            atoms = read(f"MACE-MP-0/{model}/final.xyz")
        elif model in ["MACE-ICE13", "MACE-ICE13-1"]:
            atoms = read(f"{model}/final.xyz")
        elif model == "n2p2":
            atoms = read("n2p2/01_relax-positions/final.pdb")
        return atoms

    atoms = atoms_reader(model.value)
    return atoms, atoms_reader


@app.cell
def __():
    # view(atoms)
    return


@app.cell
def __(atoms, mo):
    mo.md(
        f"""
        ## Geometria
        |valore|letteratura|simulato|errore|
        |------|-----------|--------|------|
        |Angolo alpha|{(alpha_ref := 5.5)}°|{(alpha_sim := atoms.get_angle(0, 3, 5).round(1))}°|{(alpha_sim - alpha_ref).round(2)}|
        |Angolo interno della molecola accettore| {(ang_int_acc_ref := 104.87)}° | {(ang_int_acc_sim := atoms.get_angle(1, 0, 2).round(2))}°|{(ang_int_acc_sim - ang_int_acc_ref).round(2)}
        |Angolo interno della molecola donore| {(ang_int_don_ref := 104.83)}° | {(ang_int_don_sim := atoms.get_angle(4, 3, 5).round(2))}°|{(ang_int_don_sim - ang_int_don_ref).round(2)}
        |Distanza O-O|291.2 pm = {(r_oo_ref := 2.912)} Å|{(r_oo_sim := atoms.get_distance(0, 3).round(2))} Å|{(r_oo_sim - r_oo_ref).round(2)}
        |Angolo beta|{(beta_ref := 124.4)}°|{(beta_sim := atoms.get_angle(1, 0, 3).round(1))}°|{(beta_sim - beta_ref).round(1)}
        """
    )
    return (
        alpha_ref,
        alpha_sim,
        ang_int_acc_ref,
        ang_int_acc_sim,
        ang_int_don_ref,
        ang_int_don_sim,
        beta_ref,
        beta_sim,
        r_oo_ref,
        r_oo_sim,
    )


@app.cell
def __(mo):
    mo.md("## Studio degli errori")
    return


@app.cell
def __(
    alpha_ref,
    ang_int_acc_ref,
    ang_int_don_ref,
    atoms_reader,
    beta_ref,
    mo,
    models,
    pd,
    r_oo_ref,
):
    _df = pd.DataFrame()

    for _model in models:
        _atoms = atoms_reader(_model)
        _new_df = pd.DataFrame(
            {
                "modello": _model,
                "alpha": abs(_atoms.get_angle(0, 3, 5) - alpha_ref).round(1),
                "ang_int_acc": abs(_atoms.get_angle(1, 0, 2) - ang_int_acc_ref).round(
                    2
                ),
                "ang_int_don": abs(_atoms.get_angle(4, 3, 5) - ang_int_don_ref).round(
                    2
                ),
                "r_oo": abs(_atoms.get_distance(0, 3) - r_oo_ref).round(2),
                "beta": abs(_atoms.get_angle(1, 0, 3) - beta_ref).round(1),
            },
            index=["modello"],
        )
        _df = pd.concat([_df, _new_df], axis=0, ignore_index=True)

    # _df.style.highlight_min(axis=0, color="lightgreen")
    # _df.style.background_gradient(axis=0, cmap="Greens_r")
    mo.ui.table(_df)
    return


if __name__ == "__main__":
    app.run()
