import marimo

__generated_with = "0.7.17"
app = marimo.App(width="medium")


@app.cell
def __():
    import matplotlib.pyplot as plt
    import numpy as np
    return np, plt


@app.cell
def __(np):
    freq1 = np.loadtxt("FREQ1")
    freq2 = np.loadtxt("FREQ2")
    freq3 = np.loadtxt("FREQ3")
    return freq1, freq2, freq3


@app.cell
def __(freq1, plt):
    plt.plot(freq1)
    plt.ylim(0,1)
    plt.gcf()
    return


if __name__ == "__main__":
    app.run()
