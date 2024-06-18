import marimo

__generated_with = "0.6.19"
app = marimo.App(width="medium")


@app.cell
def __():
    from ase import Atoms
    from ase.io import read
    import marimo as mo
    return Atoms, mo, read


@app.cell
def __(read):
    atoms = read("ICE13/Ih/POSCAR")
    return atoms,


@app.cell
def __(atoms):
    atoms.get_cell()
    return


@app.cell
def __(atoms):
    atoms.get_positions()
    return


@app.cell
def __(atoms):
    bandpath = atoms.cell.bandpath()
    return bandpath,


@app.cell
def __(bandpath):
    bandpath
    return


@app.cell
def __(bandpath):
    bandpath.special_points
    return


@app.cell
def __(bandpath, mo):
    mo.mpl.interactive(bandpath.plot())
    return


@app.cell
def __(mo, read):
    def plot_brillouin_zone(polymorph):
        atoms = read(f"ICE13/{polymorph}/POSCAR")
        bandpath = atoms.cell.bandpath()
        return mo.mpl.interactive(bandpath.plot())
    return plot_brillouin_zone,


@app.cell
def __(plot_brillouin_zone):
    plot_brillouin_zone("II")
    return


@app.cell
def __(plot_brillouin_zone):
    plot_brillouin_zone("III")
    return


@app.cell
def __(plot_brillouin_zone):
    plot_brillouin_zone("IV")
    return


@app.cell
def __(plot_brillouin_zone):
    plot_brillouin_zone("VI")
    return


@app.cell
def __(plot_brillouin_zone):
    plot_brillouin_zone("XI")
    return


@app.cell
def __(plot_brillouin_zone):
    plot_brillouin_zone("XV")
    return


@app.cell
def __(plot_brillouin_zone):
    plot_brillouin_zone("IX")
    return


@app.cell
def __(plot_brillouin_zone):
    plot_brillouin_zone("XVII")
    return


if __name__ == "__main__":
    app.run()
