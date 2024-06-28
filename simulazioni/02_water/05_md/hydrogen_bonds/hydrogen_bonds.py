import marimo

__generated_with = "0.5.1"
app = marimo.App()


@app.cell
def __():
    import pickle
    import numpy as np
    np.set_printoptions(linewidth=100)
    import pandas as pd

    import matplotlib.pyplot as plt

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import wa
    return mda, np, pd, pickle, plt, wa


if __name__ == "__main__":
    app.run()
