import marimo

__generated_with = "0.3.2"
app = marimo.App()


@app.cell
def __():
    import phonopy
    return phonopy,


@app.cell
def __(phonopy):
    ph = phonopy.load("POSCAR-auto-band-medium-d3-mace_mp.yaml")
    return ph,


if __name__ == "__main__":
    app.run()
