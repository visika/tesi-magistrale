import marimo

__generated_with = "0.2.13"
app = marimo.App()


@app.cell
def __(mo):
    mo.md(
        """
        # Binding energy

        Il modello è MACE-MP-0 medium con dispersione.

        L'altro è n2p2.
        """
    )
    return


@app.cell
def __():
    import pandas as pd
    return pd,


@app.cell
def __(pd):
    df = pd.read_csv('MACE-MP-0/energies.csv')
    # Sort the dataframe according to distance
    df = df.sort_values(by='distance')
    df
    return df,


@app.cell
def __(df):
    df.plot(x="distance", y="energy")
    return


@app.cell
def __():
    # Read value from file
    f = open(
        "../../01_molecule/01.8_molecule_converge_fmax_parallel_dispersion/medium/1e-8/e_gas.txt",
        "r",
    )
    e_gas_ev = float(f.read())
    # e_gas_ev = -14.170072284496875
    e_gas_ev
    return e_gas_ev, f


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
    df.plot(x="distance", y="binding_energy", ylim=(-0.4, 0.1), xlim=(2, 7))
    return


@app.cell
def __(df):
    df["binding_energy_kjmol"] = df["binding_energy"] * 96.4916
    df.plot(
        x="distance",
        y="binding_energy_kjmol",
        ylim=(-35, 0),
        xlabel="$r_{OO}$",
        xlim=(2, 7),
    )
    return


@app.cell
def __(pd):
    df_n2p2 = pd.read_csv('n2p2/energies.csv')
    df_n2p2 = df_n2p2.sort_values(by='distance')
    df_n2p2.plot(x="distance", y="energy")
    e_gas_n2p2_ev = -468.46641416393885
    df_n2p2['binding_energy'] = df_n2p2['energy'] - 2 * e_gas_n2p2_ev
    df_n2p2["binding_energy_kjmol"] = df_n2p2["binding_energy"] * 96.4916
    df_n2p2.plot(
        x="distance",
        y="binding_energy_kjmol",
        xlabel="$r_{OO}$",
        xlim=(2, 7),
    )
    return df_n2p2, e_gas_n2p2_ev


if __name__ == "__main__":
    app.run()
