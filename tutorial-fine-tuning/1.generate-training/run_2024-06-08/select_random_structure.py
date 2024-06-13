import ase
import random
from ase.io import read, write
from ase.io.trajectory import Trajectory

# trj=read('md.xyz', index=':')
trajectory = Trajectory("md.traj", "r")

number_of_structures = 50
interval_to_sample = range(2000, len(trajectory))
random_numbers = random.sample(interval_to_sample, k=number_of_structures)
print(f"Random numbers: {random_numbers}")

random_geometries = [trajectory[j] for j in random_numbers]

write(images=random_geometries, format="extxyz", filename="data_for_train.extxyz")
