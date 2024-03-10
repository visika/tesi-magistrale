# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 13-12-2023
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html

from ase.io import read,write
import matplotlib.pyplot as plt
import numpy as np


fig,ax = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0)
BIGGER_SIZE = 16


def min_loc(a,i,l):

    d = a.get_distances(i,l,mic=True)
    return np.argmin(d), np.min(d)

dt = 5.0/1000.0

rawf="bondsangles_raw.npz"
try:
    print(f"try to read raw data from {rawf}")
    raw = np.load(rawf)
    t = raw['t']
    ang = raw['aMgCO2']
    dO1 = raw['dMgCO2']
    print("done")
except:
  print(f"no raw data in {rawf}, generate")
  frames = []
  for i in range(23,24):
    frames += read(f"mg-mof-74-co2-t{i}-nvt-600.0-traj.xyz",index="200::10")
  Mgs = [a.index for a in frames[0] if a.symbol == 'Mg']
  O1=163
  O2=164
  C =162

  t=[]
  dO1=[]
  ang = []
  for i,f in enumerate(frames):
    t += [dt*i]
    i1,min1 = min_loc(f,O1,Mgs)
    i2,min2 = min_loc(f,O2,Mgs)
    if min1 < min2:
       Mg = Mgs[i1]
       O = O1
    else:
       Mg = Mgs[i2]
       O = O2
    dO1 += [ min(min1,min2) ]
    ang += [ f.get_angle(Mg,O,C,mic=True)]
  print(f"save raw data in {rawf}")
  np.savez_compressed(rawf,t=t,aMgCO2=ang,dMgCO2=dO1)

ax[0].plot(t,dO1)
ax[0].hlines(2.23,t[0],t[-1],color='red')
ax[0].fill_between(t,2.23-0.11,2.23+0.11,color='red',alpha=0.25)

ax[0].set_ylabel("d$_{Mg-O_{CO2}}$ [$\AA$]",fontsize=BIGGER_SIZE)
ax[0].set_xlabel("time [ps]",fontsize=BIGGER_SIZE)
ax[0].set_xlim(t[0],t[-1])
ax[1].plot(t,ang, "green")
ax[1].hlines(118.62,t[0],t[-1],colors='red')
ax[1].fill_between(t,118.62-10.62,118.62+10.62,color='red',alpha=0.25)
ax[1].set_ylabel("$\measuredangle_{Mg-O_{CO2}-C_{CO2}}$ [$\degree$]",fontsize=BIGGER_SIZE)
ax[1].set_xlabel("time [ps]",fontsize=BIGGER_SIZE)
ax[1].set_xlim(t[0],t[-1])
ax[1]
ax[1].tick_params(labelsize=BIGGER_SIZE)
ax[0].tick_params(labelsize=BIGGER_SIZE)

fig.align_labels()

fig.savefig("Mg-O-C.pdf")
print(f"average Mg-O_CO2 {np.average(dO1[0:])}")
print(f"average Mg-O-C {np.average(ang[0:])}")
