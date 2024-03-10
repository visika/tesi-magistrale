# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 13-12-2023
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html

from ase.io import read,write
import matplotlib.pyplot as plt
import numpy as np


metals = ['Mg','Mn','Fe','Co','Ni','Cu','Zn']
mofs = sorted( ['4b5a943', '88883f2', '92d0491', 'f22081b', '73fcf88', '38eac89', 'c4528f3', '966fade', '807cae8',
                '2d6e25d', '6f636df', 'db02aaa', '104e204', '1cf2337'])
fig,ax = plt.subplots(len(mofs), 1, layout="constrained", sharex=True, figsize=(4,len(mofs)*3))
BIGGER_SIZE = 16


def min_loc(a,i,l):

    d = a.get_distances(i,l,mic=True)
    return np.argmin(d), np.min(d)

dt = 1000.0/1000.0

def get_bond_angles(m,mo):

  rawf=f"data/bondsangles-{m}-{mo}_raw.npz"
  try:
    print(f"try to read raw data from {rawf}")
    raw = np.load(rawf)
    t = raw['t']
    ang = raw['aMgCO2']
    dO1 = raw['dMgCO2']
    print("done")
  except:
    print(f"no raw data in {rawf}, generate")
    frames = read(f"{m}/qmof-{mo}-nvt-400.0-traj.xyz",index=":")
    Met = [a.index for a in frames[0] if a.symbol == m]
    m = len(frames[0])
    O1 = m-1
    O2 = m-2
    C = m-3

    t=[]
    dO1=[]
    ang = []
    for i,f in enumerate(frames):
      t += [dt*i]
      i1,min1 = min_loc(f,O1,Met)
      i2,min2 = min_loc(f,O2,Met)
      if min1 < min2:
        Mg = Met[i1]
        O = O1
      else:
        Mg = Met[i2]
        O = O2
      dO1 += [ min(min1,min2) ]
      ang += [ f.get_angle(Mg,O,C,mic=True)]
    print(f"save raw data in {rawf}")
    np.savez_compressed(rawf,t=t,aMgCO2=ang,dMgCO2=dO1)
  return t, dO1,ang

for i,mo in enumerate(mofs):
    for m in metals:
      try:
        t, dO1, _  = get_bond_angles(m,mo)
      except:
        continue
      ax[i].plot(t,dO1,label=f"{m}")
#ax[1].plot(t,ang, "green")
#ax[1].set_ylabel("$\measuredangle_{Mg-O_{CO2}-C_{CO2}}$ [$\degree$]",fontsize=BIGGER_SIZE)
#ax[1].set_xlabel("time [ps]",fontsize=BIGGER_SIZE)
#ax[1].set_xlim(t[0],t[-1])
#ax[1].tick_params(labelsize=BIGGER_SIZE)

    ax[i].tick_params(labelsize=BIGGER_SIZE)
    ax[i].set_ylabel(f"{mo}",fontsize=BIGGER_SIZE)
    ax[i].set_xlim(0,250)
    #ax[i].legend(loc='right', bbox_to_anchor=(1, 1))
    ax[i].legend()

ax[len(mofs)-1].set_xlabel("time [ps]",fontsize=BIGGER_SIZE)
fig.align_labels()
fig.savefig(f"bonds.pdf")
