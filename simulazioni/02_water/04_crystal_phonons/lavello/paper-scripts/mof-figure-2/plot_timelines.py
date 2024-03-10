# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 16-12-2023
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html

import numpy as np
import matplotlib.pyplot as plt



metals = ['Mg','Mn','Fe','Co','Ni','Cu','Zn']
mofs = sorted( ['4b5a943', '88883f2', '92d0491', '384ba2d', 'f22081b', '73fcf88', '38eac89', 'c4528f3', '966fade', '807cae8',
                '2d6e25d', '6f636df', 'db02aaa', '391b266', '104e204', '1cf2337'])
fig,ax = plt.subplots(len(metals), 1, sharex=True)
fig.subplots_adjust(hspace=0)

for i,m in enumerate(metals):
  for mo in mofs:
    f=f"{m}/qmof-{mo}-nvt-400.0-md.log"
    ax[i].set_ylabel(f"{m} T [K]")
    try:
      d = np.loadtxt(f,skiprows=1)
    except:
      ax[i].plot()
      continue
    ax[i].plot(d[:,1]/1000,d[:,4],label=f"{m}-{mo}")
    ax[i].set_xlabel("t [ps]")
fig.savefig("T_stability.pdf")
plt.show()
