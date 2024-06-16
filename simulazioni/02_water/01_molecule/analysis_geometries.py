import marimo

__generated_with = "0.6.19"
app = marimo.App(width="medium")


@app.cell
def __():
    import marimo as mo
    from ase.io import read
    from ase.visualize import view
    from ase import Atoms
    import pandas as pd
    return Atoms, mo, pd, read, view


@app.cell
def __(read):
    atoms = read("MACE-ICE13-1/final.xyz")
    atoms
    return atoms,


@app.cell
def __(atoms):
    atoms.get_angle(0, 2, 1)
    return


@app.cell
def __():
    file_paths = [
        "01.8_molecule_converge_fmax_parallel_dispersion/small/1e-8/final.xyz",
        "01.8_molecule_converge_fmax_parallel_dispersion/medium/1e-8/final.xyz",
        "01.8_molecule_converge_fmax_parallel_dispersion/large/1e-8/final.xyz",
        "MACE-ICE13-1/final.xyz",
    ]

    models = ["small", "medium", "large", "MACE-ICE13-1"]
    return file_paths, models


@app.cell
def __(read):
    aa = read("01.8_molecule_converge_fmax_parallel_dispersion/large/1e-8/final.xyz")
    aa[aa.symbols == "H"][0]
    return aa,


@app.cell
def __(file_paths, mo, models, pd, read):
    angles = []
    bond_lengths = []

    for _f, _m in zip(file_paths, models):
        _atoms = read(_f)
        _h_indexes = [atom.index for atom in _atoms if atom.symbol == "H"]
        _o_indexes = [atom.index for atom in _atoms if atom.symbol == "O"]

        angles.append(
            {
                "model": _m,
                "angle": _atoms.get_angle(
                    _h_indexes[0], _o_indexes[0], _h_indexes[1]
                ),
            }
        )

        bond_lengths.append(
            {
                "model": _m,
                "bond_length": _atoms.get_distance(_h_indexes[0], _o_indexes[0]),
            }
        )

    # @eisenbergWaterMolecule2005
    angles.append({"model": "Reference", "angle": 104.523})

    df = pd.DataFrame(angles)

    # Compute discrepancy
    df["discrepancy"] = df["angle"] - 104.523

    # Export to CSV
    df.to_csv("angles.csv", index=False, float_format="%.3f")

    mo.ui.table(df)
    return angles, bond_lengths, df


@app.cell
def __(bond_lengths, mo, pd):
    bond_lengths_reference = {
        "model": "Reference",
        "bond_length": 0.95718,
    }
    df_bond_lengths = pd.concat(
        [
            pd.DataFrame(bond_lengths),
            pd.DataFrame(bond_lengths_reference, index=[0]),
        ]
    ).reset_index(drop=True)

    df_bond_lengths["discrepancy"] = df_bond_lengths["bond_length"] - 0.95718

    df_bond_lengths.to_csv("bond_lengths.csv", index=False, float_format="%.5f")

    mo.ui.table(df_bond_lengths)
    return bond_lengths_reference, df_bond_lengths


if __name__ == "__main__":
    app.run()
