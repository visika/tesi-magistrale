#import "@preview/whalogen:0.2.0": ce
#import "@preview/glossarium:0.4.1": make-glossary, print-glossary, gls, glspl
#import "@preview/hydra:0.4.0": hydra
#show: make-glossary

// For glossarium links
#show link: set text(fill: blue.darken(60%))

#let title = "Properties of molecular crystals using machine learning potentials"
#set document(title: title)
#set text(font: "New Computer Modern", size: 12pt)

// Justify paragraphs but not code blocks
#set par(justify: false, first-line-indent: 1em)
#show raw.where(block: true): set par(justify: false)

#set page(
  paper: "a4",
  margin: (right: 3cm, left: 3.5cm, top: 4.5cm, bottom: 3.5cm),
  numbering: "i",
  footer: [],
)

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
    columns: 2,
    align: (left, right),
    column-gutter: 8cm,
    row-gutter: 2.5mm,
    [*Relatori*], [*Candidato*],
    [Prof. Dario Alfè], [Mariano Mollo],
    [Prof. Andrea Zen], [Matr. N94000618],
  )

  #v(5.5mm)

  #text(1.2em, "Anno Accademico 2023/2024")

  #pagebreak() ]

#pagebreak()

#set page(
  paper: "a4",
  margin: (right: 3cm, left: 3.5cm, top: 3.5cm, bottom: 3.5cm),
  numbering: "i",
  footer: none,
)

#par(justify: true)[
  = Abstract

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

#pagebreak()

#outline()

#outline(
  title: [List of Figures],
  target: figure.where(kind: image),
)

#outline(
  title: [List of Tables],
  target: figure.where(kind: table),
)

#outline(
  title: [List of Code Blocks],
  target: figure.where(kind: raw),
)

#heading(numbering: none, outlined: false, [Glossary])
#print-glossary(
  (
    (
      key: "ase",
      short: "ASE",
      long: "Atomic Simulation Environment",
    ),
    (key: "md", short: "MD", long: "Molecular Dynamics"),
    (key: "dft", short: "DFT", long: "Density Functional Theory"),
    (key: "vdw", short: "vdW", long: "van der Waals"),
    (key: "mae", short: "MAE", long: "Mean Absolute Error"),
    (key: "mlp", short: "MLP", long: "Machine Learning Potential"),
    (key: "pi", short: "PI", long: "Path Integral"),
    (key: "pbe", short: "PBE", long: "Perdew-Burke-Ernzerhof"),
    (key: "dmc", short: "DMC", long: "Diffusion Monte Carlo"),
    (key: "paw", short: "PAW", long: "Projector-Augmented plane Wave"),
    (key: "zpe", short: "ZPE", long: "Zero Point Energy"),
    (key: "tnqn", short: "T&QN", long: "Thermal and Quantum Nuclear"),
    (key: "xdm", short: "XDM", long: "Exchange-hole Dipole Moment"),
    (key: "gga", short: "GGA", long: "Generalized Gradient Approximation"),
  ),
  show-all: true,
)

#pagebreak()

// https://github.com/typst/typst/discussions/3122
// New first level headings will appear on the right page
// and there will be absolutely nothing on the blank pages

#let find-labels(name) = {
  return query(name).map(label => label.location().page())
}

#let page-header = context {
  let empty-pages = find-labels(<empty-page>)
  let new-chapters = find-labels(<new-chapter>)
  if new-chapters.len() > 0 {
    if new-chapters.contains(here().page()) [
      // _a new chapter starts on this page_
      #return
    ]

    // get the index of the next <new-chapter> label
    let new-chapter-index = new-chapters.position(page => page > here().page())
    if new-chapter-index != none {
      let empty-page = empty-pages.at(new-chapter-index)
      if empty-page < here().page() [
        // _this is an empty page to make the next chapter start on an odd page_
        #return
      ]
    }
  }

  context {
    if calc.odd(here().page()) {
      align(right, hydra(1))
    } else {
      align(left, hydra(2))
    }
    // line(length: 100%)
  }
}

#let page-footer = context {
  // since the page breaks in chapter-heading() are inserted after the <empty-page> label,
  // the selector has to look "before" the current page to find the relevant label
  let empty-page-labels = query(selector(<empty-page>).before(here()))
  if empty-page-labels.len() > 0 {
    let empty-page = empty-page-labels.last().location().page()
    // look back at the most recent <new-chapter> label
    let new-chapter = query(
      selector(<new-chapter>).before(here()),
    ).last().location().page()
    // check that there is no <new-chapter> label on the current page
    if (new-chapter != here().page()) and (empty-page + 1 == here().page()) [
      // _this is an empty page where the page number should be omitted_
      #return
    ]
  }

  let page-display = counter(page).display(here().page-numbering())
  h(1fr) + page-display + h(1fr)
}

#show heading.where(level: 1): it => [
  #[] <empty-page>
  #pagebreak(to: "even", weak: true)
  #[] <new-chapter>
  #pagebreak(to: "odd", weak: true)
  // #set text(font: "Neo Euler")
  #grid(
    columns: 2,
    inset: (x: 5mm, y: 2mm),
    align: horizon + center,
    grid.vline(position: end),
    [
      #set text(size: 48pt)
      #counter(heading).display()
    ],
    upper(it.body),
  )
  #v(2em)
]

#show heading.where(level: 2): it => upper(it)

#show outline.entry.where(level: 1): it => {
  // reverse the results of the label queries to find the last <empty-page> label for the targeted page
  // the method array.position() will always return the first one...
  let empty-pages = find-labels(<empty-page>).rev()
  let new-chapters = query(<new-chapter>).rev()
  let empty-page-index = empty-pages.position(page => page == int(it.page.text))
  let new-chapter = new-chapters.at(empty-page-index)
  link(
    new-chapter.location(),
  )[#it.body #box(width: 1fr)[#it.fill] #new-chapter.location().page()]
}

#set page(header: page-header, footer: page-footer, numbering: "1")
#counter(page).update(1)

// Qui iniziano i contenuti della tesi
#set heading(numbering: "1.1")
#set par(justify: true)
#set math.equation(numbering: "(1)")

= Introduction

#box(
  stroke: 2pt + red,
  inset: 1mm,
  [
    // TODO Riprendi introduzione dall'articolo introduttivo di MACE
    *From MACE*:
    Machine-learned force fields have transformed the atomistic modelling of materials by enabling simulations of _ab initio_ quality on unprecedented time and length scales.
    However, they are currently limited by:
    1. the significant computational and human effort that must go into development and validation of potentials for each particular system of interest; and
    2. a general lack of transferability from one chemical system to the next.
    Here, using the state-of-the-art MACE architecture we introduce a single general-purpose ML model, trained on a public database of 150k inorganic crystals, that is capable of running stable molecular dynamics on molecules and materials.
    We demonstrate the power of the MACE-MP-0 model --- and its qualitative and at times quantitative accuracy --- on a diverse set problems in the physical sciences, including the properties of solids, liquids, gases, and chemical reactions.
    The model can be applied out of the box and as a starting or "foundation model" for any atomistic system of interest and is thus a step towards democatising the revolution of ML force fields by lowering the barriers to entry.
    @batatiaFoundationModelAtomistic2023
  ],
)

// TODO Descrivi cristalli molecolari e compara l'acqua come esempio caratterizzante

Qual è il problema e cosa si va a fare.

Molecular crystals are a class of materials of great technological importance.

L'acqua, i legami etc etc... Quantum Monte Carlo, DFT, quali sono i problemi.

== Motivation

DFT potentials are hard to make and it is hard to account for all the contributions.
They consume much time and resources.

The main quantity to consider to assess the stability of a crystal is its lattice energy, $E_"latt"$, which is the energy per molecule gained upon assuming the crystal form with respect to the gas state.
It can be computed as:
$
  E_"latt" = E_"crys" - E_"gas",
$ <eq-zen_2018_1>
with $E_"crys"$ as the energy per molecule in the crystal state and $E_"gas"$ as the energy of the isolated molecule.

Typically, the computation of $E_"latt"$ is performed at zero temperature and considering only the electronic contribution, i.e., quantum nuclear effects are neglected. @beranPredictingMolecularCrystal2016
The lattice energy is not directly assessable experimentally, but it can be indirectly obtained from experimental measures of the sublimation enthalpy, $Delta_"sub" H(T)$, at a given temperature $T$, by including a (theoretically evaluated) energy contribution, $Delta_"T&QN" (T)$, accounting for thermal and quantum nuclear effects:
$
  Delta_"sub" H(T) = -E_"latt" + Delta_"T&QN" (T).
$
The evaluation of $Delta_"T&QN" (T)$ can be challenging, especially for large molecules where anharmonic contributions are important. @reillyUnderstandingRoleVibrations2013
Since both $Delta_"sub" H(T)$ and $Delta_"T&QN" (T)$ are affected by errors, accurate theoretical evaluations of $E_"latt"$ are of help for comparison.

In order to derive $Delta_"T&QN"$, we need to start from the definition of the sublimation enthalpy, $Delta_"sub" H(T)$, that is the difference between the enthalpy of the gas, $H^g (T)$, and of the crystal solid, $H^s (T)$, both at temperature $T$.
By separating the electronic ($"el"$), translational ($"trans"$), rotational ($"rot"$) and vibrational ($"vib"$) contributions; noticing that in the crystal there are no trans-rotational contributions; considering as negligible the pressure times volume term, $p V$; we have that
$
  Delta_"sub" H = E_"el"^g + E_"trans"^g + E_"rot"^g + E_"vib"^g + p V - (
    E_"el"^s + E_"vib"^s
  ),
$ <eq-zen_si_14>
where the superscript stands either for gas ($g$) or solid ($s$), and the temperature dependance has been dropped for the seek of brevity.
By assuming that the rigid rotor and ideal gas approximations are reliable (that is typically the case in the analyzed molecular systems), we have that $E_"trans"^g = 3/2 R T$, $E_"rot"^g = 3/2 R T$ if the molecule is non-linear, $E_"rot"^g = R T$ otherwise, and $p V = R T$.
Thus, @eq-zen_si_14 simplifies to
$
  & Delta_"sub" H(T) = Delta E_"el" + Delta E_"vib" (
    T
  ) + 4 R T quad & "for non-linear molecules", \
  & Delta_"sub" H(T) = Delta E_"el" + Delta_"vib" (
    T
  ) + 7 / 2 R T quad & "for linear molecules",
$
where the term $Delta E_"vib" (T)$ contains both the thermal and the quantum nuclear contributions.
Notice from <eq-zen_2018_1> that $Delta E_"el" = E_"el"^g - E_"el"^s$ is precisely the opposite of the lattice energy $E_"latt"$; thus:
$
  & Delta_"T&QN" (T) = Delta E_"vib" (
    T
  ) + 4 R T quad & "for non-linear molecules", \
  & Delta_"T&QN" (T) = Delta E_"vib" (
    T
  ) + 7 / 2 R T quad & "for linear molecules".
$ <eq-zen_si_16>

Vibrations in the solid molecular crystals can usually be separated into intra-molecular and inter-molecular vibrations, $E_"vib"^s = E_"vib"^(s,"intra") + E_"vib"^(s,"inter")$, and the stiffest intra-molecular modes are decoupled from the intermolecular modes.
Intra-molecular vibrations have similar modes and frequencies than the gas-phase molecule;
thus we can conveniently write:
$
  Delta E_"vib" = Delta E_"vib"^"relax" - E_"vib"^(s,"inter")
$
where $Delta E_"vib"^"relax" := E_"vib"^g - E_"vib"^(s,"intra")$ is the change in (intra-molecular) vibrational energy given when molecules are packed in the crystal form.

At this point, a first approach can be to do a drastic approximation (which is anyway often found in the literature), that is to assume that intra-molecular frequencies in the solid are exactly the same as in the gas phase (i.e. $Delta E_"vib"^"relax" approx 0$), then to take the high temperature limit for the inter-molecular vibrations (i.e. $E_"vib"^(s,"intra") approx 6 R T$) and to neglect any zero-point motion, yielding $Delta E_"vib" (T) approx - 6 R T$ (Dulong-Petit law).
In non-linear molecules, that would imply that $Delta_"T&QN" approx -2 R T$, that is $4.96 "kJ/mol"$ at room temperature, $T=298.15K$, and zero at $T=0K$. This is a poor approximation. @otero-de-la-rozaBenchmarkNoncovalentInteractions2012
This approximation is particularly bad for water ice. @zenFastAccurateQuantum2018[SI §12.B] @whalleyEnergiesPhasesIce1984

A more reliable approach is to calculate the vibrational energies for the solid and gas phase in the harmonic limit, considering for each frequency $omega$ a contribution
$
  epsilon(omega,T) = (planck.reduce omega) / 2 + (planck.reduce omega) / (
  exp((planck.reduce omega) / (k_B T)) - 1
  ),
$
where the first term in the right hand side accounts for the @zpe contribution and the second for the thermal one.

This yields
$
  E_"vib"^g (T) = sum_i epsilon(omega_i, T), quad E_"vib"^s (
    T
  ) = integral epsilon(omega, T) g(omega) dif omega,
$ <eq-zen_si_19>
where $omega_i$ are the frequencies of the isolated molecule,
which are $3M-6$ ($M$ is the number of atoms in the molecule) for a non-linear molecule,
and $3M-5$ for a linear molecule;
$g(omega)$ is the phonon density of states in the solid. @otero-de-la-rozaBenchmarkNoncovalentInteractions2012 @reillyUnderstandingRoleVibrations2013

Notice that whenever we employ @eq-zen_si_19 to evaluate $Delta_"T&QN"$ in @eq-zen_si_16, we are subject to errors not only coming from the harmonic approximation, but also from the limitations of the computational approaches (typically @dft) used for the evaluations of the frequencies and the phonon spectrum.
Different choices of the exchange-correlation functional in @dft can lead to differences in terms of $Delta_("T&QN")$ quite larger than $1 "kJ/mol"$.
In particular, inaccuracies on the evaluation of high frequency modes mostly affects the @zpe contribution, while low frequencies modes affect mostly the thermal contribution.

=== Dispersion interactions

Dispersion interactions are essential in a collection of active research fields in solid-state physics and chemistry, including molecular crystal packing, crystal structure prediction, surface adsorption and reactivity, and supramolecular chemistry.
The representation of dispersion interactions in @dft is not possible within local or semilocal functionals because dispersion arises from non-local correlation effects involving distant fragments in the crystal. @otero-de-la-rozaBenchmarkNoncovalentInteractions2012

The @xdm model describes the dispersion energy of two neutral fragments as the electrostatic interaction of the dipoles formed by electrons and their associated exchange holes.
The dispersion energy is added to the @dft energy
$
  E = E_"DFT" + E_"disp",
$ <eq-otero_2012_1>
where $E_"disp"$ contains the usual $R^(-6)$ leading term as well as two additional higher order atomic-pairwise terms
$
  E_"disp" = -1 / 2 sum_(i j) sum_(n = 6,8,10) (C_(n,i j)) / (
  R^n_("vdw", i j) + R^n_(i j)
  ).
$ <eq-otero_2012_2>
The fundamental objects in this equation are the inter-atomic interaction coefficients $C_(n, i j)$ that in the @xdm model are calculated exclusively from first-principles quantities using second-order perturbation theory.

All the objects above are parameter-free, except for the damping expression in <eq-otero_2012_2>.
The interatomic van der Waals radii ($R_("vdw", i j)$) control the distance at which the pairwise dispersion interactions are switched off.

Because the dispersion coefﬁcients are calculated rather than ﬁtted,
@eq-otero_2012_1 works under the assumption that the @dft functional presents a completely dispersionless behavior.
This requirement is not met by most @gga functionals,
which are sometimes too repulsive and sometimes spuriously binding,
depending on the reduced-density-gradient tail behavior of the exchange enhancement factors.

==== DFT-D3
// TODO DFT-D sempre da otero
The total DFT-D3 energy is given by @grimmeConsistentAccurateInitio2010
$
  E_"DFT-D3" = E_"KS-DFT" - E_"disp",
$
where $E_"KS-DFT"$ is the usual self-consistent KS energy as obtained from the chosen DF and $E_"disp"$ is the dispersion correction as a sum of two- and three-body energies,
$
  E_"disp" = E^((2)) + E^((3)),
$
with the most important two-body term given by
$
  E^((
    2
  )) = sum_(A B) sum_(n = 6,8,10,dots) s_n (C_n^(A B)) / (r_(A B)^n) f_(d,n) (
    r_(A B)
  ).
$

Here, the first sum is over all atom pairs in the system, $C_n^(A B)$ denotes the averaged (isotropic) $n$th-order dispersion coefficient (orders $n=6,8,10,dots$) for atom pair $A B$, and $r_(A B)$ is their internuclear distance. $f_(d,n)$ are damping functions explicitly chosen by original authors to make the model numerically stable.

=== The role of anharmonic contributions
// TODO
@rossiAnharmonicQuantumFluctuations2016

Assessing predictions of lattice energies requires careful consideration of vibrational, many-body dispersion and exact-exchange contributions. @reillyUnderstandingRoleVibrations2013

== Research question

#set quote(block: true)

Taking into account the disruptive innovation brought in the molecular dynamics field by machine learning force fields, this thesis will try to answer the following question:
#quote[what is the performance of machine learning force fields in the simulation of water?]
The following chapters will introduce the theoretical foundations and the tools needed to pursue this question.

= Theory

== DFT
== Phonons
== Neural nerworks

= Tools
Ibisco, MACE @Batatia2022mace @Batatia2022Design, ASE.

== ASE

Introducing @ase.
The `Atoms` object contains the positions of the atoms and the properties of the cell.

== MACE

= Results

Seconda parte: simulazione per qualche sistema standard in cui l'approccio analitico funziona bene

Terza parte: simulazione per qualche sistema in cui l'approccio analitico non funziona bene

== Water molecule <sec-molecule>

The first task is the optimization of the geometry of the water molecule.

=== Convergence of vibrations with respect to fmax

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
    columns: 3,
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small fmax=1e-1.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small fmax=1e-4.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small fmax=1e-8.svg"),
  ),
  caption: "Convergence of the small model with respect to fmax",
)

#figure(
  grid(
    columns: 3,
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium fmax=1e-1.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium fmax=1e-4.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium fmax=1e-8.svg"),
  ),
  caption: "Convergence of the medium model with respect to fmax",
)

#figure(
  grid(
    columns: 3,
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large fmax=1e-1.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large fmax=1e-4.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large fmax=1e-8.svg"),
  ),
  caption: "Convergence of the large model with respect to fmax",
)

==== With dispersion

Inclusion of dispersion contributions leads to minimal differences in the
converged configurations.

#figure(
  grid(
    columns: 3,
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small D fmax=1e-1.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small D fmax=1e-4.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small D fmax=1e-8.svg"),
  ),
  caption: "Convergence of the small model with respect to fmax",
)

#figure(
  grid(
    columns: 3,
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium D fmax=1e-1.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium D fmax=1e-4.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium D fmax=1e-8.svg"),
  ),
  caption: "Convergence of the medium model with respect to fmax",
)

#figure(
  grid(
    columns: 3,
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large D fmax=1e-1.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large D fmax=1e-4.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large D fmax=1e-8.svg"),
  ),
  caption: "Convergence of the large model with respect to fmax",
)

==== MACE-ICE13-1

The MACE-ICE13-1 model is fine-tuned from the MACE-MP-0 medium model with
dispersion.

Convergence of results for the MACE-ICE13-1 model is achieved for displacements below $10^(-4) angstrom$.
The imaginary frequency observable in @fig-monomer-vibrations-mace-ice13-1 corresponds to negligible energy, as can be observed in the representative output in @code-monomer-vibrations-mace-ice13-1-output and therefore does not pose a issue.

#grid(
  columns: 2,
  [#figure(
      image("simulazioni/02_water/01_molecule/Grafici/MACE-ICE13-1 fmax=1e-8.svg"),
      caption: "Frequencies calculated with the MACE-ICE13-1 model at an appropriate fmax for convergence.",
    ) <fig-monomer-vibrations-mace-ice13-1>],
  [
    #figure(
      raw(
        read("simulazioni/02_water/01_molecule/MACE-ICE13-1/H2O_delta=1e-06_summary.txt"),
      ),
      caption: [
        H2O_delta=1e-06_summary.txt
      ],
    ) <code-monomer-vibrations-mace-ice13-1-output>
  ],
)

=== Zero-point vibrational energy

The computation with MACE-ICE13-1 results in a zero-point energy of $0.565 "eV"$.
This value is compared with the reference value of $+0.575 "eV"$. @eisenbergWaterMolecule2005

=== Geometry optimization

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

#figure(image("simulazioni/02_water/01_molecule/MACE-ICE13-1/ase_gui.png"))

The value found for the H-O-H angle is 104.0°. The value found for the O-H
distance is 0.970 $angstrom$.

== Water dimer

=== Vibrations analysis

The water dimer was analyzed using different models:

#figure(
  grid(
    columns: 2,
    image("simulazioni/02_water/02_dimer/01_optimize/Grafici/small.png"),
    image("simulazioni/02_water/02_dimer/01_optimize/Grafici/medium.png"),

    image("simulazioni/02_water/02_dimer/01_optimize/Grafici/large.png"),
    image("simulazioni/02_water/02_dimer/01_optimize/Grafici/MACE-ICE13-1.png"),
  ),
  caption: [
    Structure optimization of the water dimer and vibration analysis.
  ],
)

=== Geometry optimization

The optimization was performed starting from an initial configuration
approximately reproducing the geometry represented in @dimer-structure. The BFGS
optimizer was used with the MACE-MP-0 models and the MACE-ICE13-1 model, with a
force threshold of `1e-8`. The code describing the initial geometry is available
in @init.xyz. The results comparing the optimized geometry with literature
reference are available in @dimer-geometry-table and @dimer-geometry-errors.

#figure(
  image("thesis/imgs/klopper-fig1.gif"),
  caption: [
    The equilibrium structure of the water dimer. (Image taken from
    @klopperComputationalDeterminationEquilibrium2000)
  ],
) <dimer-structure>

#figure(
  ```
  6
  # CELL(abcABC):  200.00000   200.00000   200.00000    90.00000    90.00000    90.00000  Step:          67  Bead:       0 positions{angstrom}  cell{atomic_unit}
         O  1.36346e-01 -9.67442e-01  2.40661e-01
         H -7.79751e-01 -8.69284e-01  3.06661e-02
         H  5.34721e-01 -1.44347e+00 -4.71854e-01
         O  1.42482e+00  1.38927e+00  9.26440e-01
         H  1.72033e+00  1.29319e+00  1.81419e+00
         H  9.82564e-01  5.67597e-01  7.04541e-01
  ```,
  caption: [
    init.xyz
  ],
) <init.xyz>

#let dimer_geometry_table = csv("simulazioni/02_water/02_dimer/01_optimize/geometria.csv")
#figure(
  table(
    columns: dimer_geometry_table.first().len(),
    table.header(
      [Model],
      $alpha$,
      $theta_a$,
      $theta_d$,
      $r_"OO" (angstrom)$,
      $beta$,
    ),
    ..dimer_geometry_table.flatten(),
  ),
  caption: [
    Geometry values after optimization of the dimer, obtained from different
    calculator models, and from reference
    @klopperComputationalDeterminationEquilibrium2000.
    `small`, `medium`, `large` refer to MACE-MP-0.
  ],
) <dimer-geometry-table>

#let dimer_geometry_errors = csv("simulazioni/02_water/02_dimer/01_optimize/errori.csv")
#figure(
  table(
    columns: dimer_geometry_errors.first().len(),
    table.header(
      [Model],
      $alpha$,
      $theta_a$,
      $theta_d$,
      $r_"OO" (angstrom)$,
      $beta$,
    ),
    ..dimer_geometry_errors.flatten(),
  ),
  caption: [
    Errors of the geometry parameters after optimization of the dimer, obtained from
    different calculator models, with respect to reference values
    @klopperComputationalDeterminationEquilibrium2000
    @barnettBornOppenheimerMoleculardynamicsSimulations1993
    @sprikInitioMolecularDynamics1996. `small`, `medium`, `large` refer to
    MACE-MP-0.
  ],
) <dimer-geometry-errors>

=== Binding energy

The binding energy is calculated as

$ Delta E_2 := E_2 - 2 E_1 $

where $E_2$ is the energy of the dimer, $E_1$ is the energy of one molecule.

To determine the equilibrium geometry of the dimer from the study of the binding
energy, $Delta E_2$ must be minimized with respect to the internal coordinates.
To reduce the complexity of the task, the optimization is constrained to the
single degree of freedom of the OO distance.
@klopperComputationalDeterminationEquilibrium2000

Typical values for the binding energy obtained from DFT simulations at
equilibrium position range between -20 and -12 kJ/mol (approximately -0.2 to
-0.12 eV). @sprikInitioMolecularDynamics1996 @mukhopadhyayWaterDimerII2018

#figure(
  image("simulazioni/02_water/02_dimer/02_binding_energy/binding_energy.png"),
  caption: [
    Dimer binding energy for different calculators. The vertical dashed line
    corresponds to the reference value for the equilibrium OO distance,
    corresponding to the experimental value of $2.98 angstrom$
    @dykeStructureWaterDimer1977 @mukhopadhyayWaterDimerII2018. The horizontal
    dashed line corresponds to $-0.24 "eV"$
    @curtissStudiesMolecularAssociation1979.
  ],
)

== Lattice energy

The work in this section is based on @dellapiaDMCICE13AmbientHigh2022.

#figure(
  image("thesis/imgs/dellapia2022_f1.jpeg", width: 100%),
  caption: [
    Crystalline structures of the systems contained in the DMC-ICE13 dataset.
    Image obtained from @dellapiaDMCICE13AmbientHigh2022 @salzmannAdvancesExperimentalExploration2019.
  ],
)

=== Absolute lattice energy

The physical quantity usually considered to establish the stability of a crystal
is its absolute lattice energy, which is the energy per molecule gained upon
assuming the crystal form with respect to the gas phase.
It can be computed as @dellapiaDMCICE13AmbientHigh2022
$ E_"lattice" := E_"crystal" - E_"gas", $
where $E_"crystal"$ is the energy per molecule in the crystal phase, and $E_"gas"$ is the energy of the isolated molecule.

$ E_"crystal" := E / N_(#ce("H_2O")) $

The quantity $E_"gas"$ is calculated in the same manner as in @sec-molecule,
however, with the distinction that optimization was not performed;
to align the computed results with the reference paper, the same fixed gas phase molecule was used,
referenced as Patridge 1997 from the authors.

#figure(
  image("simulazioni/02_water/03_ICE13_lattice_energies/absolute_lattice_energy.svg"),
  caption: [
    Performance of MACE for the 13 ice polymorphs considered on the absolute lattice energy, compared with reference models. @dellapiaDMCICE13AmbientHigh2022
  ],
)

#figure(
  grid(
    columns: 2,
    image("simulazioni/02_water/03_ICE13_lattice_energies/scatterplot_mace-mp-0_vs_pbe.svg"),
    image("simulazioni/02_water/03_ICE13_lattice_energies/scatterplot_mace-ice13-1_vs_revpbed3.svg"),
  ),
  caption: [
    Scatter plot comparison of MACE models versus their respective reference models.
  ],
)

=== Relative lattice energy

In @dellapiaDMCICE13AmbientHigh2022 it is noted that there is interest in capturing the relative stability of the ice polymorphs, i.e., the stability with respect to a fixed crystalline phase instead of the gas state.
This property is more relevant in, e.g., the computation of the water phase diagram.
Therefore, we assess the relative stability of the crystalline phases by computing the relative lattice energy with respect to hexagonal ice Ih.

For a general polymorph $x$,
the relative lattice energy is computed as
$
  Delta E_"lattice"^x := E_"lattice"^x - E_"lattice"^"Ih"
$
and is independent of the configuration of the monomer in the gas phase.

#figure(
  image("simulazioni/02_water/03_ICE13_lattice_energies/relative_lattice_energy.svg"),
  caption: [
    Performance of MACE for the 13 ice polymorphs considered on the relative lattice energy, compared with reference models. @dellapiaDMCICE13AmbientHigh2022
  ],
)

== Crystal phonons

=== Band structure

#image("simulazioni/02_water/04_crystal_phonons/phonopy/mace_ice13_1_s2vss3_band_structure_zoom.svg")

=== Phonons DOS

#image("simulazioni/02_water/04_crystal_phonons/phonopy/mace_ice13_1_s3_dos.svg")

=== Heat capacity

#image("simulazioni/02_water/04_crystal_phonons/phonopy/heat_capacity_all_temps.svg")

== RDF

@md

#image("simulazioni/02_water/05_md/Grafici/rdf_oo_mace-ice13-1_100ps_nbins=40.svg")

= Conclusions

= Acknowledgments
// This work has been funded by project code PIR01_00011 “IBISCo”, PON 2014-2020, for all three entities (INFN, UNINA and CNR).

#show heading.where(level: 1): it => [
  #[] <empty-page>
  #pagebreak(to: "even", weak: true)
  #[] <new-chapter>
  #pagebreak(to: "odd", weak: true)
  // #set text(font: "Neo Euler")
  // #grid(
  //   columns: 2,
  //   inset: (x: 5mm, y: 2mm),
  //   align: horizon + center,
  //   grid.vline(position: end),
  //   [
  //     #set text(size: 48pt)
  //     #counter(heading).display()
  //   ],
  // upper(it.body),
  // )
  #upper(it.body)
  #v(2em)
]

#bibliography(("Tesi magistrale.bib", "bibliography.bib"))

