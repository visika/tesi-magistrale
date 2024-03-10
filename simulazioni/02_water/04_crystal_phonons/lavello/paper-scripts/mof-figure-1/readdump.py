# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 15-12-2023
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html

import numpy as np

import sys

r = []
for i in range(1,24):
    try:
      m = open(f"t{i}/c_co2_0.dump","r")
      nxt = False
      for line in m.readlines():
        if nxt:
          r += [[float(x) for x in line.split() ]]
          nxt = False
        if line.startswith("ITEM: ATOMS"):
          nxt = True
    except:
      pass

c = np.array(r)
np.savez_compressed("density_qmlff_raw.npz",C=c)



