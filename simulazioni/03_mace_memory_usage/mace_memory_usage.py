import marimo

__generated_with = "0.8.22"
app = marimo.App(width="medium")


@app.cell
def __():
    from ase.io import read
    from ase.visualize import view
    return read, view


@app.cell
def __(read, view):
    _atoms = read(
        "/home/mariano/Progetti/tesi-magistrale/strutture/ICE13/Ih/POSCAR"
    )
    view(_atoms * 3)
    return


@app.cell
def __(read):
    # Sistema molecola d'acqua, 3 atomi: /home/mariano/Progetti/tesi-magistrale/strutture/water/POSCAR
    # Sistema dimero, 6 atomi: /home/mariano/Progetti/tesi-magistrale/strutture/dimer/init.xyz
    # Sistema cella cristallina ghiaccio Ih, 12 molecole: /home/mariano/Progetti/tesi-magistrale/strutture/ICE13/Ih/POSCAR
    # Sistema con 128 molecole d'acqua: /home/mariano/Progetti/tesi-magistrale/strutture/128_molecules/h2o-128.pdb

    # atoms = read("/home/mariano/Progetti/tesi-magistrale/strutture/water/POSCAR")
    # atoms = read("/home/mariano/Progetti/tesi-magistrale/strutture/dimer/init.xyz")
    atoms = read(
        "/home/mariano/Progetti/tesi-magistrale/strutture/ICE13/Ih/POSCAR"
    )
    # atoms = read("/home/mariano/Progetti/tesi-magistrale/strutture/128_molecules/h2o-128.pdb")

    atoms = atoms * (3, 4, 4)
    # view(atoms)
    return (atoms,)


@app.cell
def __(atoms):
    print("Numero di molecole:", len(atoms)/3)
    return


@app.cell
def __():
    import matplotlib.pyplot as plt
    return (plt,)


@app.cell
def __():
    import polars as pl
    return (pl,)


@app.cell
def __(df):
    from sklearn.linear_model import LinearRegression
    import numpy as np

    X = np.array(df["molecules"]).reshape(-1, 1)
    y = np.array(df["max_memory"]).reshape(-1, 1)

    reg = LinearRegression().fit(X, y)
    return LinearRegression, X, np, reg, y


@app.cell
def __(X, plt, reg, y):
    # Plot the data and the fit together
    plt.plot(X, reg.predict(X), color="red", zorder=1)
    plt.scatter(X, y, zorder=2)
    return


@app.cell
def __(mo):
    mo.md(
        r"""
        # Simulazioni su Ibisco con la GPU

        Usa il comando:

        ```sh
        ./usr/bin/time --format="%U user + %S kernel - %M KB max memory" python script.py
        ```
        """
    )
    return


@app.cell
def __():
    import marimo as mo
    return (mo,)


@app.cell
def __(pl):
    df_ibisco = pl.DataFrame(
        {
            "molecules": [1.0, 2, 12, 128, 324, 96, 256],
            "user + kernel time": [
                6.47 + 7.69,
                6.39 + 8.01,
                6.21 + 7.65,
                7.21 + 8.57,
                7.77 + 8.66,
                6.97 + 8.14,
                7.18 + 8.42,
            ],
            "max_memory_KB": [
                1162832,
                1149524,
                1167988,
                1333676,
                1534948,
                1314772,
                1458112,
            ],
        }
    )
    return (df_ibisco,)


@app.cell
def __(pl):
    # Import data from file
    df_cuda = pl.read_csv(
        "risultati.txt",
        separator=" ",
        has_header=False,
        new_columns=["molecules", "user + kernel time", "max_memory_KB"],
    )
    return (df_cuda,)


@app.cell
def __(KB_to_GB, df_cuda, df_ibisco, pl):
    df_cuda_sanitized = df_cuda.with_columns(
        pl.col("user + kernel time")
        .map_elements(eval, return_dtype=pl.datatypes.Float64)
        .alias("user + kernel time")
    )

    # Merge df_cuda_sanitized and df_ibisco
    df_merged = df_ibisco.vstack(df_cuda_sanitized).with_columns(
        pl.col("max_memory_KB").map_batches(KB_to_GB).alias("max_memory_GB")
    )
    return df_cuda_sanitized, df_merged


@app.cell
def __():
    def KB_to_GB(KB):
        return KB / 1024 / 1024
    return (KB_to_GB,)


@app.cell
def __(mo):
    mo.md(r"""# Importazione risultati su CPU""")
    return


@app.cell
def __(KB_to_GB, pl):
    df_cpu = pl.read_csv(
        "../04_mace_usage_cpu/risultati.txt",
        separator=" ",
        has_header=False,
        new_columns=["molecules", "user + kernel time", "max_memory_KB"],
    )

    df_cpu_sanitized = df_cpu.with_columns(
        pl.col("user + kernel time")
        .map_elements(eval, return_dtype=pl.datatypes.Float64)
        .alias("user + kernel time")
    ).with_columns(
        pl.col("max_memory_KB").map_batches(KB_to_GB).alias("max_memory_GB")
    )
    return df_cpu, df_cpu_sanitized


@app.cell
def __():
    from bokeh.plotting import figure, show
    from bokeh.io import export_svg
    return export_svg, figure, show


@app.cell
def __(df_cpu_sanitized, df_merged, figure):
    _p = figure(
        title="Execution time on CPU and GPU",
        x_axis_label="Molecules",
        y_axis_label="user+kernel time (s)",
        # y_axis_type="log",
    )
    # p.y_range.start = 0
    _p.scatter(
        df_cpu_sanitized["molecules"],
        df_cpu_sanitized["user + kernel time"],
        color="blue",
        legend_label="CPU",
        marker="circle",
    )
    _p.scatter(
        df_merged["molecules"],
        df_merged["user + kernel time"],
        color="red",
        legend_label="GPU",
        marker="triangle",
    )
    _p.legend.location = "top_left"
    _p.legend.title = "Device"
    # export_svg(_p, filename="execution_time_cpu_gpu.svg")
    _p
    return


@app.cell
def __(df_cpu_sanitized, df_merged, figure):
    _p = figure(
        title="Max memory usage on CPU and GPU",
        x_axis_label="Molecules",
        y_axis_label="Max memory (GB)",
        # y_axis_type="log",
    )
    # p.y_range.start = 0
    _p.scatter(
        df_cpu_sanitized["molecules"],
        df_cpu_sanitized["max_memory_GB"],
        color="blue",
        legend_label="CPU",
        marker="circle",
    )
    _p.scatter(
        df_merged["molecules"],
        df_merged["max_memory_GB"],
        color="red",
        legend_label="GPU",
        marker="triangle",
    )
    _p.legend.location = "top_left"
    _p.legend.title = "Device"
    _p
    return


@app.cell
def __(df_merged, figure):
    _p = figure(
        title="Execution time on GPU",
        x_axis_label="Molecules",
        y_axis_label="user+kernel time (s)",
        # y_axis_type="log",
    )
    _p.y_range.start = 0
    _p.scatter(
        df_merged["molecules"],
        df_merged["user + kernel time"],
        color="red",
        legend_label="GPU",
        marker="triangle",
    )
    _p.legend.location = "bottom_right"
    _p.legend.title = "Device"
    _p
    return


@app.cell
def __(df_merged, figure):
    _p = figure(
        title="Max memory usage on GPU",
        x_axis_label="Molecules",
        y_axis_label="Max memory (GB)",
        # y_axis_type="log",
    )
    _p.y_range.start = 0
    _p.scatter(
        df_merged["molecules"],
        df_merged["max_memory_GB"],
        color="red",
        legend_label="GPU",
        marker="triangle",
    )
    _p.legend.location = "bottom_right"
    _p.legend.title = "Device"
    _p
    return


if __name__ == "__main__":
    app.run()
