import marimo

__generated_with = "0.2.5"
app = marimo.App()


@app.cell
def __():
    import pandas as pd
    return pd,


@app.cell
def __(pd):
    df_1 = pd.read_csv('energies.csv')
    df_2 = pd.read_csv('energies_2.csv')
    df = pd.concat([df_1, df_2])
    # Sort the dataframe according to distance
    df = df.sort_values(by='distance')
    df
    return df, df_1, df_2


@app.cell
def __(df):
    df.plot(x="distance", y="energy")
    return


@app.cell
def __():
    e_gas_ev = -14.170072284496875
    return e_gas_ev,


@app.cell(hide_code=True)
def __():
    import marimo as mo
    mo.md(
        """
        La binding energy si calcola come

        \[ \Delta E_2 = E_2 - 2 E_1 \]

        dove \( E_2 \) è l'energia del dimero ed \( E_1 \) è l'energia della singola molecola.
        """
    )
    return mo,


@app.cell
def __(df, e_gas_ev):
    # Calculate the binding energy as a new column
    df['binding_energy'] = df['energy'] - 2 * e_gas_ev
    df.plot(x="distance", y="binding_energy", ylim=(-0.4, 0.1))
    return


@app.cell
def __(df):
    df["binding_energy_kjmol"] = df["binding_energy"] * 96.4916
    df.plot(x="distance", y="binding_energy_kjmol", ylim=(-35, 0), xlabel="$r_{OO}$")
    return


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()