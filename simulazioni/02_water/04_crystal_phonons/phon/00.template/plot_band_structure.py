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
def __(pd):
    freq1 = pd.read_csv(
        "FREQ1", sep="[ ]+", header=None, engine="python", index_col=0
    )
    freq2 = pd.read_csv(
        "FREQ2", sep="[ ]+", header=None, engine="python", index_col=0
    )
    freq3 = pd.read_csv(
        "FREQ3", sep="[ ]+", header=None, engine="python", index_col=0
    )
    return freq1, freq2, freq3


@app.cell
def __(freq1, mo, plt):
    plt.plot(freq1, color="green")
    plt.ylim(-0.2, 5)
    # plt.savefig("freq1_s=3.svg")
    mo.mpl.interactive(plt.gcf())
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
