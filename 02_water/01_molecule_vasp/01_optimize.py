#!/usr/bin/env python3
from ase.io import Trajectory, write
from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS

calc = Vasp(command="mpirun /lustre/home/tccourse/vasp46-da/vasp", xc="pbe")

atoms = Atoms("H2O", positions=[(0, 0, 0), (0, 1, 0), (0, 0, 1)], pbc=True)
atoms.center(vacuum=3.0)
atoms.calc = calc

opt = BFGS(atoms)
traj = Trajectory("optimization.traj", "w", atoms)
opt.attach(traj)
opt.run(fmax=0.0001)
write("final.pdb", atoms)

potentialenergy = atoms.get_potential_energy()

f = open("e_gas.txt", "w")
f.write(str(potentialenergy))
f.close()
