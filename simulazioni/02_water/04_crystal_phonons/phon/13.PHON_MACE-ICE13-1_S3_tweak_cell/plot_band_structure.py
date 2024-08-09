import marimo

__generated_with = "0.7.19"
app = marimo.App(width="medium")


@app.cell
def __():
    import matplotlib.pyplot as plt
    import pandas as pd
    import marimo as mo
    return mo, pd, plt


@app.cell
def __(mo):
    directory = mo.ui.file_browser(selection_mode="directory", multiple=False, restrict_navigation=True)
    return directory,


@app.cell
def __(directory):
    directory
    return


@app.cell
def __(freq1, labels_gupta, mo, plt, ticks):
    plt.plot(freq1, color="green")
    plt.ylim(-0.2, 5)
    _ax = plt.gca()
    _ax.set_xticks(ticks, labels=labels_gupta)
    plt.grid(axis="x")
    # plt.savefig("freq1_s=3.svg")
    mo.mpl.interactive(plt.gcf())
    return


@app.cell
def __(directory):
    chosen_dir = directory.value[0].name
    return chosen_dir,


@app.cell(hide_code=True)
def __(chosen_dir, pd):
    freq1 = pd.read_csv(
        chosen_dir + "/FREQ1",
        sep="[ ]+",
        header=None,
        engine="python",
        index_col=0,
    )
    freq2 = pd.read_csv(
        chosen_dir + "/FREQ2",
        sep="[ ]+",
        header=None,
        engine="python",
        index_col=0,
    )
    freq3 = pd.read_csv(
        chosen_dir + "/FREQ3",
        sep="[ ]+",
        header=None,
        engine="python",
        index_col=0,
    )
    return freq1, freq2, freq3


@app.cell(hide_code=True)
def __():
    from ase.io import read
    atoms = read("../00.template/POSCAR.Ih")
    labels_gupta = ["G", "A", "K", "H", "M", "L", "G"]
    bandpath_gupta = atoms.cell.bandpath(labels_gupta, density=0)
    path_gupta = [bandpath_gupta.special_points.get(key) for key in labels_gupta]
    return atoms, bandpath_gupta, labels_gupta, path_gupta, read


@app.cell
def __():
    ticks = [0.00000, 0.06911, 0.18009, 0.24920, 0.33082, 0.39993, 0.50206]
    return ticks,


@app.cell
def __(bandpath_gupta):
    bandpath_gupta.get_linear_kpoint_axis()
    return


@app.cell
def __(freq2, plt):
    plt.plot(freq2)
    return


@app.cell
def __(freq3, plt):
    plt.plot(freq3)
    return


@app.cell
def __(freq1, plt):
    plt.plot(freq1)
    return


if __name__ == "__main__":
    app.run()
