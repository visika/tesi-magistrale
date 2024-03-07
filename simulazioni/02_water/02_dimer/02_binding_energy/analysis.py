import marimo

__generated_with = "0.2.13"
app = marimo.App(layout_file="layouts/analysis.grid.json")


@app.cell
def __(mo):
    mo.md(
        """
        # Binding energy

        I modelli presi in considerazione sono MACE-MP-0 medium con dispersione,
        MACE-ICE13 e n2p2.

        I valori tipici del minimo della binding energy del dimero sono
        tra i -20 e i -12 kJ/mol.
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
    df['binding_energy'] = df['energy'] - 2 * e_gas_ev
    df["binding_energy_kjmol"] = df["binding_energy"] * 96.4916
    df.plot(
        x="distance",
        y="binding_energy_kjmol",
        ylim=(-35, 0),
        xlabel="$r_{OO}$",
        ylabel="Binding energy (kJ/mol)",
        xlim=(2, 7),
        title="MACE-MP-0 medium D"
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
        title="n2p2"
    )
    return df_n2p2, e_gas_n2p2_ev


@app.cell
def __(mo):
    mo.md("## Studio di MACE-ICE13")
    return


@app.cell
def __(pd):
    df_ice13 = pd.read_csv("MACE-ICE13/energies.csv")
    df_ice13 = df_ice13.sort_values(by="distance")
    df_ice13.plot(x="distance", y="energy")
    _f = open(
        "../../01_molecule/MACE-ICE13/dispersion=True/1e-8/e_gas.txt",
        "r",
    )
    e_gas_ice13_ev = float(_f.read())
    # e_gas_ice13_ev = -14.285827697209658
    df_ice13["binding_energy"] = df_ice13["energy"] - 2 * e_gas_ice13_ev
    df_ice13["binding_energy_kjmol"] = df_ice13["binding_energy"] * 96.4916
    df_ice13.plot(
        x="distance",
        y="binding_energy_kjmol",
        xlabel="$r_{OO}$",
        ylabel="Binding energy (kJ/mol)",
        xlim=(2, 7),
        ylim=(-25, 0),
        title="MACE-ICE13 D",
    )
    return df_ice13, e_gas_ice13_ev


@app.cell
def __():
    # View the geometries
    from ase.visualize import view
    from ase.io import read
    import glob
    filenames = glob.glob("MACE-ICE13/final_dist=*.xyz")
    filenames = sorted(filenames)
    filenames
    atoms = [read(f) for f in filenames]
    # view(atoms)
    return atoms, filenames, glob, read, view


if __name__ == "__main__":
    app.run()
