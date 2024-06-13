import marimo

__generated_with = "0.6.17"
app = marimo.App(width="medium")


@app.cell
def __():
    import polars as pl

    filename = "/home/mariano/Progetti/tesi-magistrale/tutorial-fine-tuning/3.fine-tune/results/MACE-mmollo-0_run-1_train.txt"
    _data = pl.read_csv(filename, has_header=False)
    _data
    return filename, pl


@app.cell
def __():
    import pandas as pd
    import json

    # Load the data
    file_path = "/home/mariano/Progetti/tesi-magistrale/tutorial-fine-tuning/3.fine-tune/results/MACE-mmollo-0_run-1_train.txt"
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Parse the JSON Objects
    _data = [json.loads(line) for line in lines]

    # Organize the data
    df = pd.DataFrame(_data)

    # Drop NaN
    df = df.dropna()

    df
    return df, file, file_path, json, lines, pd


@app.cell
def __(df):
    # Analyze the data
    df.describe()
    return


@app.cell
def __(df):
    # Plot loss over epochs
    import matplotlib.pyplot as plt
    import os

    plt.figure(figsize=(7, 7))
    plt.plot(df["epoch"], df["loss"], label="Loss")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.title("Loss over Epochs")
    plt.legend()
    # Logscale
    plt.yscale("log")

    # Remove frames
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().spines["bottom"].set_visible(False)

    # plt.xlim(left=0)

    os.makedirs("analysis", exist_ok=True)
    # plt.savefig("analysis/loss_over_epochs.svg")
    plt.gcf()
    return os, plt


@app.cell
def __(df, os, plt):
    # Plot mae_e_per_atom over epochs
    plt.figure(figsize=(7, 7))
    plt.plot(df["epoch"], df["mae_e_per_atom"], label="MAE of energy per atom")
    plt.xlabel("Epoch")
    plt.ylabel("MAE of energy per atom")
    plt.title("MAE of energy per atom over Epochs")
    plt.legend()
    # Logscale
    plt.yscale("log")

    # Remove frames
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().spines["bottom"].set_visible(False)

    # plt.xlim(left=0)

    # plt.grid()
    os.makedirs("analysis", exist_ok=True)
    # plt.savefig("analysis/mae_e_per_atom_over_epochs.svg")
    plt.gcf()
    return


if __name__ == "__main__":
    app.run()
