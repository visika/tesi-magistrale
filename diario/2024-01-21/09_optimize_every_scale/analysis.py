import marimo

__generated_with = "0.1.79"
app = marimo.App()


@app.cell
def __():
    import pandas as pd
    return pd,


@app.cell
def __(pd):
    U0 = pd.read_csv("../08_analysis/V_vs_U0.csv")
    U0["U0"] = U0["U0"] - U0["U0"].min()
    U0.plot(x="V", y="U0", marker="o", title="U0 vs V")
    return U0,


@app.cell
def __(U0, pd, plt):
    U00 = pd.read_csv("point_zero_energy.txt", header=None)
    plt.plot(U0["V"], U00)
    return U00,


@app.cell
def __():
    import glob
    filenames = glob.glob("helmholtz*")
    # Sort naturally
    import re
    def atoi(text):
        return int(text) if text.isdigit() else text
    def natural_keys(text):
        return [atoi(c) for c in re.split("(\d+)", text)]
    filenames.sort(key=natural_keys)
    filenames
    return atoi, filenames, glob, natural_keys, re


@app.cell
def __(filenames, pd, plt):
    for f in filenames:
        df_h = pd.read_csv(f, header=None, names=["T", "F"])
        plt.plot(df_h["T"], df_h["F"], label=f)
    # plt.legend()
    plt.ylim(-469,-468)
    plt.show()
    # for f in filenames:
    #     df_h = pd.read_csv(f, header=None, names=["T", "F"])
    #     plt.plot(df_h["T"], df_h["F"], label=f)
    # plt.show()
    return df_h, f


@app.cell
def __(U0, filenames, pd):
    # Per una specifica temperatura, bisogna iterare sui volumi
    U100 = []
    for i, filename in enumerate(filenames):
        df = pd.read_csv(filename, header=None, names=["T", "F"])
        h = df[df["T"] == 100]["F"].values[0]
        U100.append(U0.iloc[i]["U0"] + h)
    U100 = U100 - min(U100)
    U200 = []
    for i, filename in enumerate(filenames):
        df = pd.read_csv(filename, header=None, names=["T", "F"])
        h = df[df["T"] == 200]["F"].values[0]
        U200.append(U0.iloc[i]["U0"] + h)
    U200 = U200 - min(U200)
    U290 = []
    for i, filename in enumerate(filenames):
        df = pd.read_csv(filename, header=None, names=["T", "F"])
        h = df[df["T"] == 290]["F"].values[0]
        U290.append(U0.iloc[i]["U0"] + h)
    U290 = U290 - min(U290)
    return U100, U200, U290, df, filename, h, i


@app.cell
def __(U0, U100, U200, U290):
    # Plot
    import matplotlib.pyplot as plt
    plt.plot(U0["V"], U0["U0"], marker="o", label="T = 0")
    plt.plot(U0["V"], U100, marker="o", label="T = 100")
    plt.plot(U0["V"], U200, marker="o", label="T = 200")
    plt.plot(U0["V"], U290, marker="o", label="T = 290")
    plt.xlabel("V")
    plt.ylabel("U")
    # plt.xlim(300, 500)
    # plt.ylim(-5, 30)
    plt.legend()
    plt.show()
    return plt,


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
