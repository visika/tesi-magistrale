#!/usr/bin/env python3
import os
from ase import units
from ase.io import read
from mace.calculators.mace import MACECalculator
from ase.md.langevin import Langevin
from datetime import datetime


atoms = read(
    filename="/ibiscostorage/mmollo/tesi-magistrale/strutture/128_molecules/h2o-128.pdb"
)

calculator = MACECalculator(
    "/ibiscostorage/VirtualMatterLab/MACE-ICE13/MACE-ICE13-1.model",
    default_dtype="float64",
    device="cuda",
)

atoms.calc = calculator

path = "MD-NVT"
os.makedirs(path, exist_ok=True)

# La durata totale della simulazione è data dal timestep per il numero di
# interazioni effettuate dal metodo run(). La temperatura nell'articolo di
# Skinner è attorno a 297.15 K.
#
# Considerando un timestep di 0.5 fs, per simulare 10 ps ho bisogno di:
# 10 * 1000 fs / 0.5 fs = 20000 iterazioni.
dyn = Langevin(
    atoms=atoms,
    timestep=0.5 * units.fs,
    temperature_K=297.15,
    friction=0.01 / units.fs,
    trajectory=f"{path}/molecular_dynamics.traj",
    loginterval=10,
)

t = datetime.now()
print(f"Molecular dynamics started at {t}")

dyn.run(10000)

dt = datetime.now() - t
print(f"Molecular dynamics finished after {dt}")
