# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 13-12-2023
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html

# -*- coding: utf-8 -*-

from ase.io import read
import numpy as np


rawf="density_raw.npz"
print(f"no raw data in {rawf}, generate")
t = []
for i in range(1,24):
  t += read(f"mg-mof-74-co2-t{i}-nvt-600.0-traj.xyz",index="::2")

iC=162
_ = [f.wrap() for f in t]
C = np.vstack([ f[iC].position for f in  t])
abc=t[0].get_cell_lengths_and_angles()
print(f"save raw data in {rawf}")
np.savez_compressed("density_raw.npz",C=C,abc=abc)
