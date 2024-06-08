import ase
from ase.io import read, write
import os

images = read("data_for_train.extxyz", ":")

for i, image in enumerate(images):
    geo = image.copy()
    os.makedirs(f"{i+1}", exist_ok=True)
    write(filename=f"{i+1}/POSCAR", images=geo, format="vasp")
