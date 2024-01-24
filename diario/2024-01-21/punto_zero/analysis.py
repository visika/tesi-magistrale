import marimo

__generated_with = "0.1.79"
app = marimo.App()


@app.cell
def __():
    import pandas as pd
    return pd,


@app.cell
def __(pd):
    # Plot data in the csv file
    pd.read_csv("V_vs_U0.csv").plot(
        x="V", y="U0", marker="o", title="Potential energy (eV)"
    )
    return


if __name__ == "__main__":
    app.run()
