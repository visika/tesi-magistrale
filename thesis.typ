= Properties of molecular crystals using machine learning potentials

= Water molecule

The first task is the optimization of the geometry of the water molecule.

== Convergence of vibrations with respect to fmax

To optimize the geometry we have to choose an optimizer. In the following, BFGS was chosen:

```python
from ase.optimize import BFGS
opt = BFGS(atoms, logfile="optimization.log", trajectory="optimization.traj")
opt.run(fmax=1e-8, steps=1000)

for d in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]:
    vib = Vibrations(atoms, delta=d, nfree=4, name=f"vib_delta={d}")
    vib.run()
    vib.summary(log=f"H2O_delta={d}_summary.txt")
    vib.clean()
```

#figure(
  grid(
    columns: 3, image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small fmax=1e-1.svg",
    ), image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small fmax=1e-4.svg",
    ), image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small fmax=1e-8.svg",
    ),
  ), caption: "Convergence of the small model with respect to fmax",
)

#figure(
  grid(
    columns: 3, image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium fmax=1e-1.svg",
    ), image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium fmax=1e-4.svg",
    ), image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium fmax=1e-8.svg",
    ),
  ), caption: "Convergence of the medium model with respect to fmax",
)

#figure(
  grid(
    columns: 3, image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large fmax=1e-1.svg",
    ), image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large fmax=1e-4.svg",
    ), image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large fmax=1e-8.svg",
    ),
  ), caption: "Convergence of the large model with respect to fmax",
)

=== With dispersion

Inclusion of dispersion contributions leads to minimal differences in the
converged configurations.

#figure(
  grid(
    columns: 3, image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small D fmax=1e-1.svg",
    ), image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small D fmax=1e-4.svg",
    ), image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small D fmax=1e-8.svg",
    ),
  ), caption: "Convergence of the small model with respect to fmax",
)

#figure(
  grid(
    columns: 3, image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium D fmax=1e-1.svg",
    ), image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium D fmax=1e-4.svg",
    ), image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium D fmax=1e-8.svg",
    ),
  ), caption: "Convergence of the medium model with respect to fmax",
)

#figure(
  grid(
    columns: 3, image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large D fmax=1e-1.svg",
    ), image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large D fmax=1e-4.svg",
    ), image(
      "simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large D fmax=1e-8.svg",
    ),
  ), caption: "Convergence of the large model with respect to fmax",
)

=== MACE-ICE13-1

#grid(columns: (2fr, 1fr),
    [The MACE-ICE13-1 model is fine-tuned from the MACE-MP-0 medium model with dispersion.],
    figure(
        image(
            "simulazioni/02_water/01_molecule/Grafici/MACE-ICE13-1 fmax=1e-8.svg",
        ), caption: "Convergence of the MACE-ICE13-1 model with respect to fmax",
    )
)

== Geometry optimization

The reference value for the H-O-H angle is 104.5°. @PhysicalChemistryWater2020
The reference value for the O-H distance is 0.96 $angstrom$. @PhysicalChemistryWater2020

= Water dimer

= Acknowledgments
This work has been funded by project code PIR01_00011 “IBISCo”, PON 2014-2020, for all three entities (INFN, UNINA and CNR).
#bibliography("bibliography.yml")
