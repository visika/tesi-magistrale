#!/usr/bin/env python3
# Compute radial distribution function (RDF)

import numpy as np
import matplotlib.pyplot as plt

from ase.geometry import analysis

# I/O to save dynamics trajectory
from ase.io.trajectory import Trajectory

traj = Trajectory("../04_molecular_dynamics/MD-NVT/simulation.traj", "r")

rdf = []
rmax = 7
nbins = 80

for i in range(len(traj)):
    geo = traj[i].copy()
    # Keep only oxygens for further analysis, discard hydrogens
    # https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.numbers
    del geo[geo.numbers == 1]
    data = analysis.get_rdf(atoms=geo, rmax=rmax, nbins=nbins, no_dists=True)
    rdf.append(data)

# Average RDF
av_rdf = 0
for i in range(len(rdf)):
    av_rdf += rdf[i] / len(rdf)

# Plot RDF
exp = np.loadtxt("rdf-exp", usecols=(0, 1))
fig, ax = plt.subplots(figsize=(6, 4))

x = np.arange(0, rmax, rmax / nbins)

ax.plot(x, av_rdf, linestyle="-", lw=2.5, c="blue", label="mlp")
ax.plot(exp[:, 0], exp[:, 1], linestyle="--", lw=2.5, c="red", label="exp")

ax.tick_params(direction="out", width=1.5, length=6)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.ylabel("$g_\mathrm{OO}(r)$", fontname="arial", fontsize=14, labelpad=8)
plt.xlabel("distance [Å]", fontname="arial", fontsize=14, labelpad=8)
plt.title("Radial Distribution Function", fontname="arial", fontsize=14)

plt.grid(ls="--", alpha=0.3)
plt.legend(fontsize=14)
plt.tight_layout()

plt.savefig("rdf.svg")
