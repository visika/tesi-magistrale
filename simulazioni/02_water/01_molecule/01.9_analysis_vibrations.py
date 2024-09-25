import marimo

__generated_with = "0.8.20"
app = marimo.App(
    app_title="Studio delle vibrazioni della molecola d'acqua",
    layout_file="layouts/01.9_analysis_vibrations.grid.json",
)


@app.cell
def __():
    import marimo as mo
    import pandas as pd
    import matplotlib.pyplot as plt
    import glob
    import os
    return glob, mo, os, pd, plt


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
    return (immaginari_a_negativi,)


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
    return (read_frequencies,)


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
    return (get_summaries,)


@app.cell
def __(mo):
    models = ["small", "medium", "large", "MACE-ICE13", "MACE-ICE13-1"]
    model = mo.ui.dropdown(options=models, value="small", label="Modello")

    fmaxs = ["1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"]
    fmax = mo.ui.dropdown(options=fmaxs, value="1e-4", label="fmax")

    dispersion = mo.ui.checkbox(value=False, label="Dispersione")

    save = mo.ui.checkbox(value=False, label="Save")

    ylabel = mo.ui.checkbox(value=True, label="Y label")

    xsize = mo.ui.number(value=8, start=1, stop=20, step=1, label="X size")
    ysize = mo.ui.number(value=6, start=1, stop=20, step=1, label="Y size")


    mo.hstack([model, fmax, dispersion, ylabel, xsize, ysize, save])
    return (
        dispersion,
        fmax,
        fmaxs,
        model,
        models,
        save,
        xsize,
        ylabel,
        ysize,
    )


@app.cell
def __(mo):
    fontsize = mo.ui.number(value=15, start=1, stop=20, step=1, label="Font size")
    ticksize = mo.ui.number(start=12, stop=20, step=1, label="Ticks size")
    mo.hstack([ticksize, fontsize])
    return fontsize, ticksize


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
def __(
    df,
    dispersion,
    fmax,
    fontsize,
    mo,
    model,
    os,
    plt,
    save,
    ticksize,
    ylabel,
):
    # plt.figure(figsize=(xsize.value, ysize.value))
    plt.subplots(layout="constrained")

    # Set labels size
    plt.rcParams.update({"font.size": ticksize.value})

    groups = df.groupby("#")
    for _name, _group in groups:
        plt.plot(_group["delta"], _group["cm^-1"], marker="x", label=_name)
    plt.xscale("log")
    plt.yscale("symlog", linthresh=1e-1)
    plt.ylim(-1e4, 1e4)
    plt.legend(
        ncol=3,
        loc="lower left",
    )

    if model.value == "MACE-ICE13-1":
        title_fmax = "fmax=1e-8"
        model_string = model.value
    else:
        title_fmax = f"fmax={fmax.value}"
        model_string = f"{'MACE-MP-0' if model.value in ['small', 'medium', 'large'] else ''} {model.value}"


    dispersion_string = f"{' D' if dispersion.value else ''}"

    plt.title(
        title_string := f"H2O vib. modes "
        + model_string
        + dispersion_string
        + ", "
        f"{title_fmax}",
        fontsize=fontsize.value,
    )
    plt.xlabel("Displacement (Å)", fontsize=fontsize.value)

    if ylabel.value:
        plt.ylabel("Frequency ($\mathrm{cm}^{-1}$)")
    # plt.grid()

    # Hide the top of the frame
    plt.gca().spines[["top", "bottom", "left", "right"]].set_visible(False)

    if save.value:
        save_folder = "Grafici"
        os.makedirs(save_folder, exist_ok=True)
        save_path = (
            f"{save_folder}/{model_string}{dispersion_string} {title_fmax}.svg"
        )
        plt.savefig(save_path)
        print(f"Saved to {save_path}")

    # Show an interactive marimo plot
    mo.mpl.interactive(plt.gcf())
    return (
        dispersion_string,
        groups,
        model_string,
        save_folder,
        save_path,
        title_fmax,
        title_string,
    )


@app.cell
def __(mo):
    mo.md(rf"## Studio delle energie di punto zero")
    return


@app.cell
def __():
    import re

    def read_zero_point_energy(filename):
        with open(
            filename,
            "r",
        ) as _f:
            zpe_line = _f.readlines()[-1]
            zpe = re.findall(r"[-+]?(?:\d*\.*\d+)", zpe_line)[0]
            return zpe
    return re, read_zero_point_energy


@app.cell
def __(get_summaries, mo, pd, read_zero_point_energy):
    _fmax = "1e-8"
    _dispersion = True
    _zpes = []
    _models = ["small", "medium", "large", "MACE-ICE13-1"]

    # Le ZPE calcolate sono riguardanti le sole vibrazioni armoniche
    # Quindi bisogna usare come riferimento la ZPE della molecola di acqua in approssimazione armonica
    # La ZPE_harmonic vale 0.585 eV, dalla reference Barone2004

    zpe_harmonic_eV = 0.585

    for _model in _models:
        _summary_filename = get_summaries(_model, _fmax, _dispersion)[0]
        _zpes.append(float(read_zero_point_energy(_summary_filename)))
    _df = pd.concat(
        [
            pd.DataFrame({"model": _models, "zpe": _zpes}),
            # @eisenbergWaterMolecule2005
            pd.DataFrame({"model": "Expt.", "zpe": zpe_harmonic_eV}, index=[0]),
        ]
    )

    # Compute the errors
    _df["error"] = round(_df["zpe"] - zpe_harmonic_eV, 3)

    # Export to CSV
    _df.to_csv("zero_point_energies.csv", index=False, float_format="%.3f")

    mo.ui.table(_df)
    return (zpe_harmonic_eV,)


@app.cell
def __(mo):
    mo.md(rf"## Studio dei tre modi normali di vibrazione")
    return


@app.cell
def __(get_summaries, pd, read_frequencies):
    _models = ["small", "medium", "large", "MACE-ICE13-1"]


    def le_mie_frequenze(model):
        _fname = get_summaries(model, "1e-8", True)[0]
        _df = read_frequencies(_fname)
        return _df.tail(3)


    ni_1_reference = 3656.65
    ni_2_reference = 1594.59
    ni_3_reference = 3755.79

    # Tre discorsi separati per ciascuna frequenza
    ni_1 = pd.DataFrame({"model": "Reference", "cm^-1": ni_1_reference}, index=[0])
    ni_2 = pd.DataFrame({"model": "Reference", "cm^-1": ni_2_reference}, index=[0])
    ni_3 = pd.DataFrame({"model": "Reference", "cm^-1": ni_3_reference}, index=[0])

    for _model in _models:
        _df = le_mie_frequenze(_model)
        ni_1 = pd.concat(
            [
                ni_1,
                pd.DataFrame(
                    {"model": _model, "cm^-1": _df["cm^-1"].iloc[-2]}, index=[0]
                ),
            ],
            ignore_index=True,
        )
        ni_2 = pd.concat(
            [
                ni_2,
                pd.DataFrame(
                    {"model": _model, "cm^-1": _df["cm^-1"].iloc[-3]}, index=[0]
                ),
            ],
            ignore_index=True,
        )
        ni_3 = pd.concat(
            [
                ni_3,
                pd.DataFrame(
                    {"model": _model, "cm^-1": _df["cm^-1"].iloc[-1]}, index=[0]
                ),
            ],
            ignore_index=True,
        )

    # Calcolo delle discrepanze
    ni_1["error"] = ni_1["cm^-1"] - ni_1_reference
    ni_2["error"] = ni_2["cm^-1"] - ni_2_reference
    ni_3["error"] = ni_3["cm^-1"] - ni_3_reference

    # Esportazione in CSV
    ni_1.to_csv("ni_1.csv", index=False, float_format="%.2f")
    ni_2.to_csv("ni_2.csv", index=False, float_format="%.2f")
    ni_3.to_csv("ni_3.csv", index=False, float_format="%.2f")
    return (
        le_mie_frequenze,
        ni_1,
        ni_1_reference,
        ni_2,
        ni_2_reference,
        ni_3,
        ni_3_reference,
    )


@app.cell
def __(mo, ni_1, ni_2, ni_3):
    mo.hstack([ni_1, ni_2, ni_3])
    return


@app.cell
def __():
    from ase.vibrations import Vibrations
    from ase.io import read

    from mace.calculators import mace_mp

    calculator = mace_mp(model="medium", device="cpu", default_dtype="float64")

    h2o = read("MACE-ICE13-1/final.xyz")
    h2o.calc = calculator

    vib = Vibrations(h2o)
    vib.clean()
    vib.run()
    vib.write_mode(-1)
    vib.write_mode(-2)
    vib.write_mode(-3)
    vib.summary()
    return Vibrations, calculator, h2o, mace_mp, read, vib


@app.cell
def __(mo):
    mo.md(rf"## Studio dei valori della reference in approssimazione armonica")
    return


@app.cell
def __():
    omega_1_reference = 3832.17
    omega_2_reference = 1648.47
    omega_3_reference = 3942.53
    return omega_1_reference, omega_2_reference, omega_3_reference


@app.cell
def __(
    le_mie_frequenze,
    omega_1_reference,
    omega_2_reference,
    omega_3_reference,
    pd,
):
    # Tre discorsi separati per ciascuna frequenza
    omega_1 = pd.DataFrame(
        {"model": "Reference", "cm^-1": omega_1_reference}, index=[0]
    )
    omega_2 = pd.DataFrame(
        {"model": "Reference", "cm^-1": omega_2_reference}, index=[0]
    )
    omega_3 = pd.DataFrame(
        {"model": "Reference", "cm^-1": omega_3_reference}, index=[0]
    )

    _models = ["small", "medium", "large", "MACE-ICE13-1"]
    for _model in _models:
        _df = le_mie_frequenze(_model)
        omega_1 = pd.concat(
            [
                omega_1,
                pd.DataFrame(
                    {"model": _model, "cm^-1": _df["cm^-1"].iloc[-2]}, index=[0]
                ),
            ],
            ignore_index=True,
        )
        omega_2 = pd.concat(
            [
                omega_2,
                pd.DataFrame(
                    {"model": _model, "cm^-1": _df["cm^-1"].iloc[-3]}, index=[0]
                ),
            ],
            ignore_index=True,
        )
        omega_3 = pd.concat(
            [
                omega_3,
                pd.DataFrame(
                    {"model": _model, "cm^-1": _df["cm^-1"].iloc[-1]}, index=[0]
                ),
            ],
            ignore_index=True,
        )

    # Calcolo delle discrepanze
    omega_1["error"] = omega_1["cm^-1"] - omega_1_reference
    omega_2["error"] = omega_2["cm^-1"] - omega_2_reference
    omega_3["error"] = omega_3["cm^-1"] - omega_3_reference
    return omega_1, omega_2, omega_3


@app.cell
def __(mo, omega_1, omega_2, omega_3):
    mo.hstack([omega_1, omega_2, omega_3])
    return


@app.cell
def __(omega_1, omega_2, omega_3, os, pd):
    _df0 = omega_1["model"]
    _df1 = omega_1["cm^-1"]
    _df2 = omega_2["cm^-1"]
    _df3 = omega_3["cm^-1"]
    _df_final = pd.concat([_df0, _df1, _df2, _df3], axis=1)
    _df_final.columns = ["Model", "ω1", "ω2", "ω3"]
    _save_folder = "Analysis"
    os.makedirs(_save_folder, exist_ok=True)
    _df_final.to_csv(f"{_save_folder}/omega.csv", index=False, float_format="%.2f")
    _df_final
    return


@app.cell
def __(omega_1, omega_2, omega_3, pd):
    # Now the same for the errors
    _df0 = omega_1["model"]
    _df1 = omega_1["error"]
    _df2 = omega_2["error"]
    _df3 = omega_3["error"]
    _df_final = pd.concat([_df0, _df1, _df2, _df3], axis=1)
    _df_final.columns = ["Model", "ω1", "ω2", "ω3"]
    _save_folder = "Analysis"
    _df_final.to_csv(
        f"{_save_folder}/omega_errors.csv", index=False, float_format="%.2f"
    )
    _df_final
    return


if __name__ == "__main__":
    app.run()
