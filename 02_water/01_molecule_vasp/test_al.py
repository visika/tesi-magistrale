#!/usr/bin/env python3
from ase.io import Trajectory, write
from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS
from ase.build import bulk

al = bulk("Al")
mydir = "bandstructure"
calc = Vasp(
    command="mpirun /ibiscostorage/VirtualMatterLab/vasp_6.3.0_bin/vasp_std",
    kpts=(2, 2, 2),
    directory=mydir,
)
al.calc = calc
al.get_potential_energy()
kpts = {"path": "WGX", "npoints": 30}
calc.set(isym=0, icharg=11, kpts=kpts)
print(al.get_potential_energy())
bs = calc.band_structure()
print(bs)

# atoms = Atoms("H2O", positions=[(0, 0, 0), (0, 1, 0), (0, 0, 1)], pbc=True)
# atoms.center(vacuum=3.0)
# atoms.calc = calc

# opt = BFGS(atoms)
# traj = Trajectory("optimization.traj", "w", atoms)
# opt.attach(traj)
# opt.run(fmax=0.0001)
# write("final.pdb", atoms)

# potentialenergy = atoms.get_potential_energy()

# f = open("e_gas.txt", "w")
# f.write(str(potentialenergy))
# f.close()
