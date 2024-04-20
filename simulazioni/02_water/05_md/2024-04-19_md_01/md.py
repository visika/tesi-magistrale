#!/usr/bin/env python3

import os
from ase import units
from ase.io import read
from mace.calculators.mace import MACECalculator
from ase.md.langevin import Langevin


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

dyn = Langevin(
    atoms=atoms,
    timestep=0.5 * units.fs,
    temperature_K=300.0,
    friction=0.01 / units.fs,
    trajectory=f"{path}/molecular_dynamics.traj",
    loginterval=10,
)

dyn.run(1000)
