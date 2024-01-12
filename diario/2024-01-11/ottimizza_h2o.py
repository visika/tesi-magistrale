#!/usr/bin/env python3
# https://wiki.fysik.dtu.dk/ase/gettingstarted/tut02_h2o_structure/h2o.html
from ase import Atoms
from ase.optimize import BFGS
from gpaw import GPAW

atoms = Atoms("HOH", positions=[(0, 0, -1), (0, 1, 0), (0, 0, 1)])
atoms.center(vacuum=3.0)

calc = GPAW(mode="lcao", basis="dzp", txt="gpaw.txt")
atoms.calc = calc

opt = BFGS(atoms, trajectory="opt.traj")
opt.run(fmax=0.01)

print(atoms.get_angle(0, 1, 2))
