# lavello


collection of random but useful scripts, before they made it somewhere else.

practically we showcase [mace mlip](https://github.com/ACEsuit/mace) usage in some typical computational physics tasks on top of various python packages, ase, phonopy and
so on.

are not used for production runs, the main scope is proof of concept.

it downloads for you the model if you use small, medium or large keywords consistent with this paper: https://arxiv.org/abs/2401.00096
you can pass a model of your choice by passing the path to it.

## content

- md.py simple but complex md driver for using mace ml ff with ASE:
  - with or without Grimme d3 dispersion correction (pass --d3 to activate)
  - nve, nvt, npt
  - geometry minimization atoms and cell
  - heating
  - --help for all supported options
- eos.py simple equation of state driver, using ase
- phonons.py computes phonons, force constants, using phonopy... quick dos/pdos and bands also calculated and ploted
- neb.py simple nudge elastic bands calculations using ase


all scripts have command line arguments to control various options, use --help to see what is possible

## examples

### md.py

 - heating (starts by minimizing your coordinates and cell, then starts a classical md (npt ensemble) at T=1 K up to T=1000K in
   steps of 10 K, each temperature step is 1000 fs)

```bash
    sys="UiO66" && python3 ./md.py --barostat_time 2500 --system structs/$sys.cif --model medium  --minimize --optimize_cell --minimize_fmax 0.05  --temp_start 1.0 --temp_end 1000.0 --temp_step 10 --temp_time 1000 --ensemble npt --system_name $sys --restart_rotate
```

 - simple nvt md, using Langevin dynamics, Nose Hoover a la Melchionna is possible.

```bash
   sys="UiO66"  &&  python3 ./md.py --system structs/$sys.cif --temperature 300 --ensemble nvt --seed $RANDOM --steps 2000 --system_name $sys --traj_every 1000 --output_every 100 --timestep 1.0 --model medium --restart_rotate --device cpu
```

note cpu will be quite slow for a system of few hundred atoms use a gpu for orders of magnitude faster times.

- geometry optimisation (positions and cell)

```bash
  sys="UiO66"  &&  python3 ../md.py --system $sys.cif --system_name $sys --minimize --minimize_fmax 0.1 --optimize_cell  --model medium --device cpu
```

- single point calculation

```bash
  sys="UiO66"  &&  python3 ../md.py --system $sys.cif  --model medium  --device cpu
```

### eos example

```bash
   python3 ./eos.py --system structs/Ag.xyz --model large --device cpu
```

### phonons

 - it uses phonopy and seekpath so install them in addition to mace

```bash
   python3 ./phonons.py --model medium --system structs/Si.cif  --device cpu --grid 20 20 20
```

you will get a file *Si-ase-medium.yaml* that contains the force constants for further processing


### neb

- run a symple neb between two states this is same example as here -
  https://github.com/materialsvirtuallab/matcalc/blob/main/examples/LiFePO4-NEB.ipynb (note, mace gets barriers closes to the DFT
PBE values) - reproduces some data from this paper https://doi.org/10.1039/C5TA05062F


```bash
   python3 ./neb.py --A structs/LiFePO4_bc.cif --B structs/LiFePO4_c.cif --device cpu --model medium
```
