import marimo

__generated_with = "0.8.22"
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
    plt.savefig("freq1_s=2.svg")
    plt.ylim(-0.2, 5)
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


@app.cell
def __(pd):
    dos = pd.read_csv(
        "DOS", sep="[ ]+", header=None, engine="python", index_col=0
    )
    return (dos,)


@app.cell
def __(dos):
    dos
    return


@app.cell
def __(dos):
    dos.plot(xlim=(0, 12))
    return


@app.cell
def __(pd):
    smeared_dos = pd.read_csv(
        "Smeared_DOS_Data.csv"
    )
    return (smeared_dos,)


@app.cell
def __(smeared_dos):
    smeared_dos
    return


@app.cell
def __(dos):
    dos.iloc[:,0]
    return


@app.cell
def __(smearedd_dos):
    len(smearedd_dos)
    return


@app.cell
def __(plt, smeared_dos):
    plt.plot(smeared_dos["Energy"], smeared_dos["Smeared_DOS"])
    plt.xlim(0, 12)
    plt.gca()
    return


@app.cell
def __(dos):
    dos.index
    return


@app.cell
def __(dos, pd, smearedd_dos):
    # Create a new DataFrame with smeared values
    smeared_data = pd.DataFrame({
        'Energy': dos.index,
        'Smeared_DOS': smearedd_dos
    })
    return (smeared_data,)


@app.cell
def __(plt, smeared_data):
    plt.plot(smeared_data["Energy"], smeared_data["Smeared_DOS"])
    plt.xlim(0, 12)
    plt.grid()
    plt.gca()
    return


@app.cell
def __(dos):
    from scipy.ndimage import gaussian_filter1d

    # Apply a small Gaussian smearing to the second column (DOS values)
    sigma = 1.0 # A small Gaussian sigma value
    smearedd_dos = gaussian_filter1d(dos.iloc[:, 0], sigma=sigma)
    return gaussian_filter1d, sigma, smearedd_dos


if __name__ == "__main__":
    app.run()
