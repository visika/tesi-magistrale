from ase import io
from ase.build import add_adsorbate

from ase.build import molecule

import argparse
import pathlib
import numpy as np

cli=argparse.ArgumentParser()

cli.add_argument(
  "--system",
  type=str,
  default='',
  help = ' provide system to add co2. default %(default)s'
)

args = cli.parse_args()
system = args.system


d = system.split("/")[0]
stem = pathlib.Path(system).stem

s = io.read(system)
print(d,stem)
h = s.get_cell_lengths_and_angles()
for m in ['CO2']:
    mol=molecule(m)
    add_adsorbate(s,mol, height=-0.5*h[2],position=(0,0.65*h[1]*np.cos(np.radians(30.0))))
    io.write(f"{d}/{stem}-co2.cif",s,format="cif")


