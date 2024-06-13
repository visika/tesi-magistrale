import marimo

__generated_with = "0.6.17"
app = marimo.App(width="medium")


@app.cell
def __():
    from ase.io.trajectory import Trajectory
    return Trajectory,


@app.cell
def __(Trajectory):
    trajectory = Trajectory("md.traj", "r")
    return trajectory,


@app.cell
def __():
    from ase.visualize import view
    return view,


@app.cell
def __(trajectory, view):
    view(trajectory[-1])
    return


@app.cell
def __(trajectory):
    len(trajectory)
    return


@app.cell
def __(trajectory):
    import random
    interval = range(9990, len(trajectory))
    random_numbers = random.sample(interval, k=10)
    random_numbers
    return interval, random, random_numbers


@app.cell
def __(trajectory):
    # Get the energies from the trajectory
    energies = [atoms.get_potential_energy() for atoms in trajectory]
    # Plot the energies
    import matplotlib.pyplot as plt
    plt.plot(energies)
    plt.xlabel("Frame")
    plt.ylabel("Energy (eV)")
    plt.show()
    return energies, plt


@app.cell
def __(trajectory):
    trajectory
    return


if __name__ == "__main__":
    app.run()
