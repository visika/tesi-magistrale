from ase.io import read
from mace.calculators import mace_mp

# Sistema molecola d'acqua, 3 atomi: /ibiscostorage/mmollo/tesi-magistrale/strutture/water/POSCAR
# Sistema dimero, 6 atomi: /ibiscostorage/mmollo/tesi-magistrale/strutture/dimer/init.xyz
# Sistema cella cristallina ghiaccio Ih, 12 molecole: /ibiscostorage/mmollo/tesi-magistrale/strutture/ICE13/Ih/POSCAR
# Sistema con 128 molecole d'acqua: /ibiscostorage/mmollo/tesi-magistrale/strutture/128_molecules/h2o-128.pdb

# atoms = read("/ibiscostorage/mmollo/tesi-magistrale/strutture/water/POSCAR")
# atoms = read("/ibiscostorage/mmollo/tesi-magistrale/strutture/dimer/init.xyz")
atoms = read("/ibiscostorage/mmollo/tesi-magistrale/strutture/ICE13/Ih/POSCAR")
# atoms = read("/ibiscostorage/mmollo/tesi-magistrale/strutture/128_molecules/h2o-128.pdb")

# Eventuale supercella
atoms = atoms * (3, 3, 4)

calc = mace_mp(model="medium", dispersion=True, default_dtype="float64", device="cuda")
atoms.calc = calc

print("Total energy of the system:", atoms.get_total_energy())
print("Forces in the system:")
print(atoms.get_forces())

print("Numero di molecole:", len(atoms) / 3)
