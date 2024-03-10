# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 13-12-2023
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html

# -*- coding: utf-8 -*-

from ase.io import read
import numpy as np
import pathlib as p

metals = ['Mg','Mn','Fe','Co','Ni','Cu','Zn']
mofs = sorted( ['88883f2', '92d0491', 'f22081b', '73fcf88', '38eac89', 'c4528f3', '966fade', '807cae8',
                '2d6e25d', '6f636df', 'db02aaa', '104e204', '1cf2337'])
iC=0
for m in metals:
  for mo in mofs:
    try:
      fn = f"{m}/qmof-{mo}-co2-traj.xyz"
      t = read(fn,index=":")
    except:
      print(f"no data for {metal}")
      continue

    print(f"read {len(t)} frames for {m}-qmof-{mo} from {fn}")
    _ = [f.wrap() for f in t]
    C = np.vstack([ f[iC].position for f in  t])
    abc=t[0].cell.cellpar()
    data_file=f"{m}-{mo}.npz"

    print(f"save raw data in data/{data_file}")
    np.savez_compressed(f"data/{data_file}",C=C,abc=abc)
