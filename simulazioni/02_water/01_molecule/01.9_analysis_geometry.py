import marimo

__generated_with = "0.6.20"
app = marimo.App()


@app.cell
def __():
    import marimo as mo
    from ase.io import read
    from ase.visualize import view
    import matplotlib.pyplot as plt
    return mo, plt, read, view


@app.cell
def __():
    def get_results(atoms, angle_atoms=[0, 2, 1], bond_atoms=[0, 2]):
        results = {}
        # Get the H-O-H angle between atoms 0, 2, 1
        results["angle"] = atoms.get_angle(*angle_atoms)
        # Get the O-H bond length between atoms 0, 2
        results["bond_length"] = atoms.get_distance(*bond_atoms)
        return results
    return get_results,


@app.cell
def __(mo):
    mo.md("# MACE-ICE13-1")
    return


@app.cell
def __(get_results, read):
    with open("/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/01_molecule/MACE-ICE13-1/final.xyz", "r") as _f:
        atoms_mace_ice13_1 = read(_f, format="xyz")
    get_results(atoms_mace_ice13_1)
    return atoms_mace_ice13_1,


@app.cell
def __(mo):
    mo.md("# MACE-MP-0 + D3")
    return


@app.cell
def __(get_results, read):
    with open("/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/01_molecule/01.8_molecule_converge_fmax_parallel_dispersion/medium/1e-8/final.xyz", "r") as _f:
        atoms_mace_mp_0_medium_D = read(_f, format="xyz")
    get_results(atoms_mace_mp_0_medium_D, angle_atoms=[0, 1, 2], bond_atoms=[0, 1])
    return atoms_mace_mp_0_medium_D,


@app.cell
def __(get_results, read):
    for _model in ["small", "medium", "large"]:
        with open(
            f"/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/01_molecule/01.8_molecule_converge_fmax_parallel_dispersion/{_model}/1e-8/final.xyz",
            "r",
        ) as _f:
            _atoms = read(_f, format="xyz")
        print(
            f"{_model}:\t {get_results(_atoms, angle_atoms=[0, 1, 2], bond_atoms=[0, 1])}"
        )
    return


@app.cell
def __(get_results, plt, read):
    _angoli = {}
    for _fmax in ["1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"]:
        with open(
            f"/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/01_molecule/01.8_molecule_converge_fmax_parallel_dispersion/large/{_fmax}/final.xyz",
            "r",
        ) as _f:
            _atoms = read(_f, format="xyz")
        _angoli[_fmax] = get_results(
            _atoms, angle_atoms=[0, 1, 2], bond_atoms=[0, 1]
        )["angle"]
    print(_angoli)

    plt.plot([float(e) for e in _angoli.keys()], _angoli.values())

    # Logscale
    plt.xscale("log")

    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    plt.ylim(top=104.495)

    plt.xlabel("fmax")
    plt.ylabel("H-O-H angle")
    plt.title("MACE-MP-0 large + D3")

    # remove frames
    for _spine in plt.gca().spines.values():
        _spine.set_visible(False)

    plt.tight_layout()
    plt.gca()
    return


@app.cell
def __(mo):
    mo.md("# MACE-MP-0")
    return


@app.cell
def __(get_results, read):
    with open("/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/01_molecule/01.7_molecule_converge_fmax_parallel/medium/1e-8/final.xyz", "r") as _f:
        atoms_mace_mp_0_medium = read(_f, format="xyz")
    get_results(atoms_mace_mp_0_medium, angle_atoms=[0, 1, 2], bond_atoms=[0, 1])
    return atoms_mace_mp_0_medium,


@app.cell
def __(get_results, read):
    for _model in ["small", "medium", "large"]:
        with open(
            f"/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/01_molecule/01.7_molecule_converge_fmax_parallel/{_model}/1e-8/final.xyz",
            "r",
        ) as _f:
            _atoms = read(_f, format="xyz")
        print(
            f"{_model}:\t {get_results(_atoms, angle_atoms=[0, 1, 2], bond_atoms=[0, 1])}"
        )
    return


@app.cell
def __(get_results, plt, read):
    angoli = {}
    for fmax in ["1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"]:
        with open(
            f"/home/mariano/Progetti/tesi-magistrale/simulazioni/02_water/01_molecule/01.7_molecule_converge_fmax_parallel/large/{fmax}/final.xyz",
            "r",
        ) as _f:
            atoms = read(_f, format="xyz")
        angoli[fmax] = get_results(
            atoms, angle_atoms=[0, 1, 2], bond_atoms=[0, 1]
        )["angle"]
    print(angoli)

    plt.plot([float(e) for e in angoli.keys()], angoli.values())

    # Logscale
    plt.xscale("log")

    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)

    plt.xlabel("fmax")
    plt.ylabel("H-O-H angle")
    plt.title("MACE-MP-0 large")

    # remove frames
    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    # plt.savefig("Grafici/angle_convergence_mace_mp_0_large.svg")

    plt.gca()
    return angoli, atoms, fmax, spine


@app.cell
def __(angoli):
    last_value = angoli["1e-8"]
    diffs = {k: abs(v - last_value) for k, v in angoli.items()}
    diffs
    return diffs, last_value


@app.cell
def __(angoli, last_value):
    # Get relative differences
    rdiffs = {k: abs(v - last_value) / v for k, v in angoli.items()}
    rdiffs
    return rdiffs,


if __name__ == "__main__":
    app.run()
