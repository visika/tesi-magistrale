import marimo

__generated_with = "0.3.2"
app = marimo.App()


@app.cell
def __():
    import phonopy
    return phonopy,


@app.cell
def __(phonopy):
    ph = phonopy.load("kgrid=8/ase-medium-d3-mace_mp.yaml")
    return ph,


@app.cell
def __(ph):
    # for m in [2, 4, 8, 16]:
    for m in [32]:
        ph.run_mesh([m, m, m])
        
        ph.run_total_dos(freq_min=0, freq_max=14)
        ph.write_total_dos(filename=f"total_dos_mesh_{m}.dat")

    bplt = ph.plot_total_dos()
    bplt.xlim((0, 15))
    # bplt.ylim((0, 8))
    bplt.show()
    return bplt, m


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
