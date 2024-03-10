# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 13-12-2023
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html

# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import numpy as np

from matplotlib.ticker import MaxNLocator
from ase.units import kB
import matplotlib as mplt
import matplotlib.cm as cm
import scipy.cluster.vq as scv

import pathlib as p
import matplotlib as mpl
import pickle
import argparse

cli=argparse.ArgumentParser()
cli.add_argument(
  "--pdf",
  action="store_true",
  help = 'render pdf? default %(default)s'
)

cli.add_argument(
  "--no_centers",
  action="store_true",
  help = 'do not show centres? default %(default)s'
)

cli.add_argument(
  "--no_autoscale",
  action="store_true",
  help = 'no common scale default %(default)s'
)

SMALL_SIZE = 8
MEDIUM_SIZE = 18
BIGGER_SIZE = 28
cmap = plt.colormaps['viridis']
metals = ['Mg','Mn','Fe','Co','Ni','Cu','Zn']
mofs = sorted( ['88883f2', '92d0491', 'f22081b', '73fcf88', '38eac89', 'c4528f3', '966fade', '807cae8',
                '2d6e25d', '6f636df', 'db02aaa', '104e204', '1cf2337'])
#mofs = sorted( ['4b5a943', '88883f2', '92d0491', '384ba2d', 'f22081b', '73fcf88', '38eac89', 'c4528f3', '966fade', '807cae8'])
broken =[('Zn','4b5a943')]
nmet = len(metals)
nmof =len(mofs)
step = 0
args = cli.parse_args()
pdf = args.pdf
noc = args.no_centers
autos = args.no_autoscale


def heat_plot(a: np.ndarray,xlabels: list,ylabels: list,title:str, filename:str) -> None :

     hig, ax = plt.subplots(layout="constrained")
     im = ax.imshow(a,cmap=cmap)
     ax.set_xticks(np.arange(len(xlabels)), labels=xlabels,fontsize=SMALL_SIZE)
     ax.set_yticks(np.arange(len(ylabels)), labels=ylabels,fontsize=SMALL_SIZE)
     plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")
     for i in range(len(metals)):
       for j in range(len(mofs)):
          text = ax.text(j, i, f"{a[i, j]:.3f}",
                       ha="center", va="center", color="w")
     cbar = ax.figure.colorbar(im, location="top",ax=ax )
     cbar.ax.set_ylabel(title, rotation=0,va="bottom", fontsize=SMALL_SIZE)
     cbar.ax.tick_params(labelsize=SMALL_SIZE)
     cbar.ax.yaxis.set_label_coords(0.5,2)

     hig.savefig(filename)


# compute free energy
labels = ['density','free-energy']
def save_image(filename: str) -> None:

    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]

    for i,fig in enumerate(figs):
        p = PdfPages(f"{filename}-{labels[i]}.pdf")
        fig.savefig(p, format='pdf')
        p.close()


def get_density_free_enegy(C: np.ndarray ,nx: int, ny: int)  :
    H, xedges,yedges = np.histogram2d(C[step:,0],C[step:,1], bins = (nx, ny),density=True)
    p = np.zeros(H.shape)
    for i in range(nx):
      for j in range(ny):
        if H[i,j] > eps:
          p[i,j] = -kBT*np.log(H[i,j])
        else:
          p[i,j] = -kBT*np.log(eps)
      eref = kBT*np.log(eps)
      F=  eref + p
    return np.min(C,axis=0), np.max(C,axis=0), H, F, xedges, yedges

def insets(g,nmet,mofs):
  for i,mo in enumerate(mofs):
    g[nmet,i].imshow(plt.imread(f"qmof-{mo}.png") )
    g[nmet,i].axis('off')

def grid_plot(g,f,A,cent,X,Y,title,levels,extend_max):
  X = (X[:-1] + X[1:]) / 2
  Y = (Y[:-1] + Y[1:]) / 2
# free energy
  if not noc:
    s = g[i,j].scatter(cent[(m,mo)][0],cent[(m,mo)][1],s=3.0,c='b', marker='o')
  if (m,mo) in broken:
    p = g[i,j].contourf(X,Y,A, levels=levels, cmap=cmap,alpha=.75, antialiased=True)
  else:
    p = g[i,j].contourf(X,Y,A, levels=levels, cmap=cmap)
  c = p.get_cmap()
  if extend_max:
    col = cmap(np.linspace(0,1,nz))[0]
  else:
    col = cmap(np.linspace(0,1,nz))[-1]
  g[i,j].get_xaxis().set_visible(False)
  g[i,j].get_yaxis().set_visible(False)
  g[i,j].set_facecolor(col)

  if i == len(metals)-1:
      g[i,j].set_title(f"{mo}",y=-0.275,fontsize=BIGGER_SIZE)
  if j == 0:
    g[i,j].set_ylabel(f"{m}",fontsize=BIGGER_SIZE)
    g[i,j].get_yaxis().set_visible(True)
    g[i,j].set_yticks([])
  if i==j and i == 0:
    cb = f.colorbar(p,location="left",ax=g[:,0],shrink=0.75,**kwargs)
    cb.ax.tick_params(labelsize=BIGGER_SIZE,rotation=0)
    cb.ax.set_ylabel(title,fontsize=BIGGER_SIZE+4)
    cb.ax.yaxis.set_label_coords(-0.125,1.5)


dig, dx = plt.subplots(nrows=nmet+1, ncols=nmof,
        figsize=(nmof*2,(nmet+1)*2),layout="constrained",
        )
Fig, Fx = plt.subplots(nrows=nmet+1, ncols=nmof,
        figsize=(nmof*2,(nmet+1)*2),layout="constrained",
        )

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

temp = 400.0
kBT = kB*temp
eps = 1.0e-5

nx = 100
ny = 100
nz = 100


densmin = np.zeros((len(metals),len(mofs)))
densmax = np.zeros((len(metals),len(mofs)))
freemin = np.zeros((len(metals),len(mofs)))
freemax = np.zeros((len(metals),len(mofs)))

for i,m in enumerate(metals):
  for j,mo in enumerate(mofs):
    try:
      raw = np.load(f"data/{m}-{mo}.npz")
    except:
      continue
    C = raw['C']
    mi, ma, H, F, xedges, yedges =  get_density_free_enegy(C,nx,ny)
    densmin[i,j] = H.min()
    densmax[i,j] = H.max()
    freemin[i,j] = F.min()
    freemax[i,j] = F.max()

heat_plot(densmax,mofs,metals,"max projected density","density-heat.pdf")
heat_plot(freemin,mofs,metals,"min free energy [eV]","free-energy-heat.pdf")

Hmin = np.min(densmin)
Hmax = np.max(densmax)
Fmin = np.min(freemin)
Fmax = np.max(freemax)
kwargs = {'format': '%.2f'}

cent = pickle.load(open("data/centers.pkl","rb"))

levels = MaxNLocator(nbins=nz).tick_values(Hmin,Hmax)
levels2 = MaxNLocator(nbins=nz).tick_values(Fmin,Fmax)

for i,m in enumerate(metals):
  for j,mo in enumerate(mofs):
    try:
      raw = np.load(f"data/{m}-{mo}.npz")
    except:
      continue

# density
    C = raw['C']
    abc = raw['abc']
    mi, ma, H, F, X, Y =  get_density_free_enegy(C,nx,ny)
    print(f"{m}-{mo} ρ: min={H.min():.3f} max={H.max():.3f} F: min={F.min():.3f} max={F.max():.3f}")

    if not autos:
      levels = MaxNLocator(nbins=nz).tick_values(H.min(),H.max())
    grid_plot(dx,dig,H,cent,X,Y,"projected density",levels,True)
    grid_plot(Fx,Fig,F,cent,X,Y,"free energy [eV]",levels2,False)

insets(dx,len(metals),mofs)
insets(Fx,len(metals),mofs)

if pdf:
  save_image("mof-74")
dig.savefig("mof-74-density.png")
Fig.savefig("mof-74-free-energy.png")
