#import "@preview/whalogen:0.2.0": ce

#let title = "Properties of molecular crystals using machine learning potentials"
#set document(title: title)
#set text(font: "New Computer Modern", size: 12pt)
#set page(
  paper: "a4",
  margin: (right: 3cm, left: 3.5cm, top: 4.5cm, bottom: 3.5cm),
)
#set heading(numbering: "1.")

#[ // Title page
  #set align(center)

  #text(1.5em, [*UNIVERSITÀ DEGLI STUDI DI NAPOLI* \ *"FEDERICO II"*])

  #v(3mm)

  // University Logo
  #image("thesis/imgs/University_Federico_II_Logo.svg", width: 25%)

  #v(1cm)

  *Scuola Politecnica e delle Scienze di Base*

  *Area Didattica di Scienze Matematiche Fisiche e Naturali*

  #v(8mm)

  *Dipartimento di Fisica "Ettore Pancini"*

  #v(20mm)

  _Laurea Magistrale in Fisica_

  #v(5mm)

  #text(1.5em, title)

  #v(25mm)

  #grid(
    columns: 2, align: (left, right), column-gutter: 8cm, row-gutter: 2.5mm, [*Relatori*], [*Candidato*], [Prof. Dario Alfè], [Mariano Mollo], [Prof. Andrea Zen], [Matr. N94000618],
  )

  #v(5.5mm)

  #text(1.2em, "Anno Accademico 2023/2024")

  #pagebreak() ]

#set page(
  paper: "a4", margin: (right: 3cm, left: 3.5cm, top: 3.5cm, bottom: 3.5cm),
)

#par(
  justify: true,
)[
  *Abstract*

  Molecular crystals play an important role in the field of materials science,
  particularly in drug development, electronics, and renewable energy sectors.

  In this work we will study the properties of molecular crystals, using recently
  developed Machine Learning potentials to model their behaviour and
  characteristics. We will be primarily focusing on water as the initial subject,
  followed by a study of a selection of other molecular crystals.

  Traditional approaches often grapple with the trade-off between computational
  expense and accuracy. The application of Machine Learning potentials captures
  complex intermolecular interactions with a significantly reduced computational
  cost compared to traditional ab-initio methods.

  We will study the capabilities of trained Machine Learning potentials to
  accurately predict lattice energies, polymorphic behaviours, and response to
  external conditions like temperature and pressure. We will also study dynamic
  properties such as phonon spectra to complete the insight into the physical and
  chemical behaviours of molecular crystals.
]

#outline()

= Introduction

Qual è il problema e cosa si va a fare.

Molecular crystals are a class of materials of great technological importance.

L'acqua, i legami etc etc... Quantum Monte Carlo, DFT, quali sono i problemi.

== Prima parte: studio di letteratura
= Tools
Ibisco, MACE @Batatia2022mace @Batatia2022Design, ASE.

= Seconda parte: simulazione per qualche sistema standard in cui l'approccio analitico funziona bene
= Terza parte: simulazione per qualche sistema in cui l'approccio analitico non funziona bene

= Water molecule

The first task is the optimization of the geometry of the water molecule.

== Convergence of vibrations with respect to fmax

To optimize the geometry we have to choose an optimizer. In the following, BFGS
was chosen:

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

#grid(
  columns: (2fr, 1fr), [The MACE-ICE13-1 model is fine-tuned from the MACE-MP-0 medium model with
    dispersion.], figure(
    image("simulazioni/02_water/01_molecule/Grafici/MACE-ICE13-1 fmax=1e-8.svg"), caption: "Frequencies calculated with the MACE-ICE13-1 model at converged fmax",
  ),
)

== Geometry optimization

The reference value for the H-O-H angle is 104.5°. @PhysicalChemistryWater2020
The reference value for the O-H distance is 0.96 $angstrom$.
@PhysicalChemistryWater2020

The simulation was performed using the fine-tuned MACE-ICE13-1 model.

```python
calculator = MACECalculator(
    "/ibiscostorage/VirtualMatterLab/MACE-ICE13/MACE-ICE13-1.model",
    default_dtype="float64",
    device="cuda",
)
opt = BFGS(atoms, logfile="optimization.log", trajectory="optimization.traj")
opt.run(fmax=1e-8, steps=1000)
write("final.xyz", images=atoms)
```

The final geometry was analyzed through the visualizer integrated in ase:

```sh
ase gui final.xyz
```

#image("simulazioni/02_water/01_molecule/MACE-ICE13-1/ase_gui.png")

The value found for the H-O-H angle is 104.0°. The value found for the O-H
distance is 0.970 $angstrom$.

= Water dimer

== Vibrations analysis

The water dimer was analyzed using different models:

#figure(
  grid(
    columns: 2, image("simulazioni/02_water/02_dimer/01_optimize/Grafici/small.png"), image("simulazioni/02_water/02_dimer/01_optimize/Grafici/medium.png"), image("simulazioni/02_water/02_dimer/01_optimize/Grafici/large.png"), image("simulazioni/02_water/02_dimer/01_optimize/Grafici/MACE-ICE13-1.png"),
  ), caption: [
    Structure optimization of the water dimer and vibration analysis.
  ],
)

== Geometry parameters

#let dimer_geometry_errors = csv("simulazioni/02_water/02_dimer/01_optimize/errori.csv")
#figure(
  table(
    columns: dimer_geometry_errors.first().len(), ..dimer_geometry_errors.flatten(),
  ), caption: [
    Errors of the geometry parameters of the dimer, obtained from different
    calculator models, with respect to reference values
    @barnettBornOppenheimerMoleculardynamicsSimulations1993
    @sprikInitioMolecularDynamics1996.
  ],
)

== Binding energy

= Acknowledgments
// This work has been funded by project code PIR01_00011 “IBISCo”, PON 2014-2020, for all three entities (INFN, UNINA and CNR).
#bibliography(("bibliography.yml", "bibliography.bib"))
