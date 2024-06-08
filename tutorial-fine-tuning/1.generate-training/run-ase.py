import numpy as np
from ase.io import read
from datetime import datetime

# from mace.calculators.mace import MACECalculator
# from ase.md.langevin import Langevin  # NVT molecular dynamics
from ase.md.nptberendsen import NPTBerendsen
from ase.md import MDLogger
import ase.units as units
from mace.calculators import mace_mp

# mace=MACECalculator('/rds/user/fd385/hpc-work/work/MolCrys/ICE/it_0/models/MACE_model_large_swa.model', default_dtype='float64',device='cuda')
print("Loading calculator")
calculator = mace_mp(
    model="medium", dispersion=False, default_dtype="float32", device="cuda"
)
atoms = read("POSCAR", 0)
atoms.calc = calculator
print("Finished loading calculator")

# dyn = Langevin(
#     atoms,
#     timestep=0.5 * units.fs,
#     temperature_K=100.0,  # temperature in K
#     friction=0.01 / units.fs,
#     trajectory="md.traj",
# )

print("Setting dynamics")
dyn = NPTBerendsen(
    atoms,
    timestep=0.5 * units.fs,
    temperature_K=100.0,
    pressure_au=1.01325 * units.bar,
    # Compressibility depends on temperature
    compressibility_au=0.9471 * 1e-10 / units.Pascal,
    taut=100 * units.fs,
    taup=500 * units.fs,
    trajectory="md.traj",
)
# 0.5 ps = 500 fs

dyn.attach(
    MDLogger(dyn, atoms, "md.log", header=True, stress=False, mode="a"), interval=4
)
print("Finished setting dynamics")
t = datetime.now()
print(f"Starting dynamics at {t}")
dyn.run(10000)
print(f"Finished dynamics after {datetime.now() - t}")
