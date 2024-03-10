# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 13-12-2023
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html

# -*- coding: utf-8 -*-

from ase.io import read
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import numpy as np

from matplotlib.ticker import MaxNLocator
from ase.units import kB
import matplotlib as mplt
import matplotlib.cm as cm
import scipy.cluster.vq as scv

SMALL_SIZE = 6
MEDIUM_SIZE = 8
BIGGER_SIZE = 12
steps = 0

def colormap2arr(arr,cmap):
    # http://stackoverflow.com/questions/3720840/how-to-reverse-color-map-image-to-scalar-values/3722674#3722674
    gradient=cmap(np.linspace(0.0,1.0,100))
    # Reshape arr to something like (240*240, 4), all the 4-tuples in a long list...
    arr2=arr.reshape((arr.shape[0]*arr.shape[1],arr.shape[2]))

    # Use vector quantization to shift the values in arr2 to the nearest point in
    # the code book (gradient).
    code,dist=scv.vq(arr2,gradient)

    # code is an array of length arr2 (240*240), holding the code book index for
    # each observation. (arr2 are the "observations".)
    # Scale the values so they are from 0 to 1.
    values=code.astype('float')/gradient.shape[0]

    # Reshape values back to (240,240)
    values=values.reshape(arr.shape[0],arr.shape[1])
    values=values[::-1]
    return values


# compute free energy
labels = ['density','free-energy']
def save_image(filename):

    # PdfPages is a wrapper around pdf
    # file so there is no clash and
    # create files with no error.

    # get_fignums Return list of existing
    # figure numbers
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]

    # iterating over the numbers in list
    for i,fig in enumerate(figs):
        p = PdfPages(f"{filename}-{labels[i]}.pdf")

        # and saving the files
        fig.savefig(p, format='pdf')

    # close the object
        p.close()

def get_density_free_enegy(C,nx,ny):
    H, xedges,yedges = np.histogram2d(C[steps:,0],C[steps:,1], bins = (nx, ny),density=True)
    p = np.zeros(H.shape)
    for i in range(nx):
      for j in range(ny):
        if H[i,j] > eps:
          p[i,j] = -kBT*np.log(H[i,j])
        else:
          p[i,j] = -kBT*np.log(eps)
      eref = kBT*np.log(eps)
      F=  eref + p
    xcenters = (xedges[:-1] + xedges[1:]) / 2
    ycenters = (yedges[:-1] + yedges[1:]) / 2
    xmi, xma  = xcenters.min(), xcenters.max()
    ymi, yma = ycenters.min(), ycenters.max()
    X = xcenters - (xma+xmi)*0.5
    Y = ycenters - (yma+ymi)*0.5
    return np.min(C,axis=0), np.max(C,axis=0), H, F, X, Y

temp = 600.0
kBT = kB*temp
eps = 1.0e-5

rawf="density_raw.npz"
qrawf="density_qmlff_raw.npz"
try:
  print(f"try to read raw data from {rawf}")
  raw = np.load(rawf)
  qraw = np.load(qrawf)
  C = raw['C']
  qC = qraw['C']
  abc = raw['abc']
#read classical picture and make it a graph

  arr=plt.imread('uff-mof.png')
  uff=colormap2arr(arr,cm.jet)

  ny,nx = uff.shape
  xx =np.linspace(-5,5,nx)
  yy =np.linspace(-5,5,ny)
  cX, cY = np.meshgrid(xx, yy)

  cF = uff.copy()
  a = 0.03
  b = -0.14
  for i in range(ny):
    for j in range(nx):
      if abs(cF[i,j]-0.49) < 1.0e-5:
        cF[i,j] = 0.0
      cF[i,j] = a+cF[i,j]*(b-a)

  print(f"done")
except:
  print(f"no raw data in {rawf}, generate")

nx =  100
ny =  100
nz = 100

mi, ma, H, F, X, Y =  get_density_free_enegy(C,nx,ny)
mi, ma, qH, qF, qX, qY =  get_density_free_enegy(qC,nx,ny)

print(f"H: min={H.min()} max={H.max()}")
print(f"F: min={F.min()} max={F.max()}")
print(f"qH: min={qH.min()} max={qH.max()}")
print(f"qF: min={qF.min()} max={qF.max()}")
fig, ax = plt.subplots(
        layout='constrained',
        figsize=(3, 3),
        )

cmap = plt.colormaps['viridis']
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


levels = MaxNLocator(nbins=nz).tick_values(H.min(),H.max())

p = ax.contourf(X,Y,H, levels=levels, cmap=cmap)
cb=fig.colorbar(p)
cb.ax.set_ylabel('projected density [$\AA^{-2}$]')
ax.set_xlim(np.min(X), np.max(X)*0.99)
ax.set_ylim(np.min(Y), np.max(Y)*0.99)
ax.set_xlabel("x [$\AA$]")
ax.set_ylabel("y [$\AA$]")

fig1, ax1 = plt.subplots(1,3,
        layout='constrained',
        figsize=(9, 3))


levels2 = MaxNLocator(nbins=nz).tick_values(F.min(),F.max())
p = ax1[0].contourf(X,Y,F, levels=levels2, cmap=cmap)
ax1[0].set_xlim(np.min(X), np.max(X)*0.99)
ax1[0].set_ylim(np.min(Y), np.max(Y)*0.99)
ax1[0].get_xaxis().set_visible(False)
ax1[0].get_yaxis().set_visible(False)
ax1[0].set_title("MACE-MP0",loc="left",fontsize=BIGGER_SIZE)

c = p.get_cmap()
col = cmap(np.linspace(0,1,nz))[-1]

p = ax1[1].contourf(qY,qX,qF, levels=levels2, cmap=cmap,extend='max')
ax1[1].set_xlim(np.min(X), np.max(X)*0.99)
ax1[1].set_ylim(np.min(Y), np.max(Y)*0.99)
ax1[1].get_xaxis().set_visible(False)
ax1[1].get_yaxis().set_visible(False)
ax1[1].set_title("QMLFF",loc="left",fontsize=BIGGER_SIZE)
ax1[1].set_facecolor(col)


p = ax1[2].contourf(cY,cX,cF, levels=levels2, cmap=cmap,extend='max')
ax1[2].set_xlim(np.min(X), np.max(X)*0.99)
ax1[2].set_ylim(np.min(Y), np.max(Y)*0.99)
ax1[2].get_xaxis().set_visible(False)
ax1[2].get_yaxis().set_visible(False)
ax1[2].set_facecolor(col)
ax1[2].set_title("UFF",loc="left",fontsize=BIGGER_SIZE)

cb=fig1.colorbar(p)
cb.ax.tick_params(labelsize=BIGGER_SIZE)
cb.ax.set_ylabel('free energy [eV]',fontsize=BIGGER_SIZE)
fig1.savefig("free_energy.png")
save_image("Mg-MOF-74-CO2")
