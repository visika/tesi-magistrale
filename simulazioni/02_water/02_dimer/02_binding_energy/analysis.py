import marimo

__generated_with = "0.9.9"
app = marimo.App(width="columns", app_title="Binding energy")


@app.cell(column=0, hide_code=True)
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
    return (pd,)


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
    return (mo,)


@app.cell
def __(pd):
    def df_binding_energy(model: str, molecule_energy):
        """
        model: name of the model
        molecule_energy: energy of the molecule in eV
        """
        df = pd.read_csv(f"{model}/energies.csv")
        df = df.sort_values(by="distance")
        df["binding_energy"] = df["energy"] - 2 * molecule_energy
        df["binding_energy_kjmol"] = df["binding_energy"] * 96.4916
        return df
    return (df_binding_energy,)


@app.cell
def __(df_binding_energy):
    with open(
        "../../01_molecule/01.8_molecule_converge_fmax_parallel_dispersion/medium/1e-8/e_gas.txt",
        "r",
    ) as _f:
        e_gas_ev = float(_f.read())
        df = df_binding_energy("MACE-MP-0", e_gas_ev)
    return df, e_gas_ev


@app.cell
def __(df_binding_energy):
    e_gas_n2p2_ev = -468.46641416393885
    df_n2p2 = df_binding_energy("n2p2", e_gas_n2p2_ev)
    return df_n2p2, e_gas_n2p2_ev


@app.cell(disabled=True)
def __(df_binding_energy):
    with open(
        "../../01_molecule/MACE-ICE13/dispersion=True/1e-8/e_gas.txt", "r"
    ) as _f:
        e_gas_ice13_ev = float(_f.read())
        df_ice13 = df_binding_energy("MACE-ICE13", e_gas_ice13_ev)
    return df_ice13, e_gas_ice13_ev


@app.cell
def __(df_binding_energy):
    with open(
        "/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/01_molecule/MACE-ICE13-1/e_gas.txt",
        "r",
    ) as _f:
        e_gas_mace_ice13_1 = float(_f.read())
        df_mace_ice13_1 = df_binding_energy("MACE-ICE13-1", e_gas_mace_ice13_1)
    return df_mace_ice13_1, e_gas_mace_ice13_1


@app.cell
def __():
    # ❯ pint-convert eV kcal/mol
    # 1 electron_volt = 23.0605478306 kcal/mol


    def eV2kcalmol(x):
        return x * 23.0605478306


    def kcalmol2eV(x):
        return x / 23.0605478306
    return eV2kcalmol, kcalmol2eV


@app.cell
def __(df, df_mace_ice13_1, eV2kcalmol, kcalmol2eV, mo):
    # Plot all the data
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(layout="constrained")

    # Draw a horizontal line
    # kjmol : -22.7
    # 22.7 kilojoule / mole = 0.235268921197 eV
    ax.axhline(y=-0.235268921197, color="r", linestyle="--", label="Equilibrium exp. value")

    # Draw a vertical line
    # angstrom : 2.98
    ax.axvline(x=2.98, color="r", linestyle="--")

    ax.plot(
        df_mace_ice13_1["distance"],
        df_mace_ice13_1["binding_energy"],
        label="MACE-ICE13-1",
        marker="s",
        markersize=5,
    )

    ax.plot(
        df["distance"],
        df["binding_energy"],
        label="MACE-MP-0 medium+D",
        marker="o",
        markersize=5,
    )

    # ax.plot(
    #     df_n2p2["distance"],
    #     df_n2p2["binding_energy"],
    #     label="n2p2",
    #     marker="^",
    #     markersize=5,
    # )

    ax.set_xlabel("$r_{OO}$ (Å)")
    ax.set_ylabel("Binding energy (eV)")

    ax.set_ylim(-0.35, 0)

    ax.legend()
    plt.title("Dimer binding energy")

    # remove the left and right spines
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    secax = ax.secondary_yaxis("right", functions=(eV2kcalmol, kcalmol2eV))
    secax.spines["right"].set_visible(False)
    secax.set_ylabel("Binding energy (kcal/mol)")

    # fig.savefig("binding_energy.svg")

    # interactive plot
    mo.mpl.interactive(plt.gcf())
    return ax, fig, plt, secax


@app.cell(disabled=True)
def __():
    # View the geometries
    from ase.visualize import view
    from ase.io import read
    import glob
    filenames = glob.glob("MACE-ICE13-1/final_dist=*.xyz")
    filenames = sorted(filenames)
    filenames
    atoms = [read(f) for f in filenames]
    view(atoms)
    return atoms, filenames, glob, read, view


@app.cell(column=1)
def __(mo):
    mo.md(r"""## Confronto con DFT per la presentazione con le slide""")
    return


@app.cell
def __():
    # Valori calcolati a r_oo = 3.0 Å
    revPBE_dimer = -28.6759816
    revPBE_monomer = -14.25742076
    return revPBE_dimer, revPBE_monomer


@app.cell
def __(revPBE_dimer, revPBE_monomer):
    revPBE_binding_energy = revPBE_dimer - 2 * revPBE_monomer
    return (revPBE_binding_energy,)


@app.cell
def __(df_binding_energy, revPBE_monomer):
    revPBE_binding_energies = df_binding_energy("VASP/varia_r", revPBE_monomer)
    revPBE_binding_energies
    return (revPBE_binding_energies,)


@app.cell
def __(df, df_mace_ice13_1, mo, plt, revPBE_binding_energies):
    _fig, _ax = plt.subplots(layout="constrained")

    _ax.plot(
        revPBE_binding_energies["distance"],
        revPBE_binding_energies["binding_energy"],
        label="revPBE",
        marker="^",
        markersize=7,
        color="green",
    )

    _ax.plot(
        df_mace_ice13_1["distance"],
        df_mace_ice13_1["binding_energy"],
        label="MACE-ICE13-1",
        marker="s",
        markersize=5,
        color="blue",
    )

    _ax.plot(
        df["distance"],
        df["binding_energy"],
        label="MACE-MP-0 medium",
        marker="o",
        markersize=5,
        color="orange",
    )

    # _ax.plot(3.0, revPBE_binding_energy, marker="^", markersize=10, label="revPBE")

    _ax.set_xlabel("$r_{OO}$ (Å)")
    _ax.set_ylabel("Energia (eV)")

    _ax.set_ylim(-0.35, 0)

    _ax.legend()
    plt.title("Energia di interazione del dimero")

    # remove the left and right spines
    _ax.spines["left"].set_visible(False)
    _ax.spines["right"].set_visible(False)
    _ax.spines["bottom"].set_visible(False)

    # _fig.savefig("binding_energy_revPBE.svg")
    mo.mpl.interactive(plt.gcf())
    return


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
