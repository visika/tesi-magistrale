#import "@preview/whalogen:0.2.0": ce
#import "@preview/glossarium:0.4.1": make-glossary, print-glossary, gls, glspl
#import "@preview/hydra:0.4.0": hydra
#import "@preview/fletcher:0.4.5" as fletcher: diagram, node, edge
#import "@preview/tablem:0.1.0": tablem
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
    (key: "csp", short: "CSP", long: "Crystal Structure Prediction"),
    (key: "pes", short: "PES", long: "Potential Energy Surface"),
    (key: "qti", short: "QTI", long: "Quantum Thermodynamic Integration"),
    (key: "vasp", short: "VASP", long: "Vienna Ab initio Simulation Package"),
    (key: "xc", short: "XC", long: "Exchange–Correlation"),
    (key: "rdf", short: "RDF", long: "Radial Distribution Function"),
    (
      key: "ccsdt",
      short: "CCSD(T)",
      long: "coupled cluster theory with singlets, doublets, and perturbative triplets",
    ),
    (key: "cbs", short: "CBS", long: "complete basis set"),
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
  #set text(font: "Neo Euler")
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

#let large_figure(content, caption: none) = figure(
  box(width: 125%, content),
  caption: caption,
)

#let large_box(content) = align(
  center,
  box(
    width: 125%,
    content,
  ),
)

= Introduction

#box(
  stroke: 2pt + green,
  inset: 1mm,
  [
    The structure and properties of molecular crystals have long been of great interest,
    not just for fundamental reasons of understanding molecular aggregation but also due to their numerous applications.
    As the preferred form of active pharmaceutical ingredients for oral administration, the dissolution and morphology of drug-molecule crystals are very important for bio-availability and processing. @blagdenCrystalEngineeringActive2007
    For these reasons alone the prediction of molecular crystal structures is of the utmost importance. @priceComputationalPredictionPharmaceutical2004
    In addition, molecular crystals can also have a wide range of optical, electronic, and mechanical properties, @reddyMechanicalPropertiesMolecular2010 which in some cases can be tuned based on environmental variables @martinsTemperaturePressureInducedProton2009 or composition. @karkiImprovingMechanicalProperties2009 @reillyUnderstandingRoleVibrations2013
  ],
)

#figure(
  image("thesis/imgs/karki2009_graphical_abstact.jpg", width: 50%),
  caption: [
    Poor mechanical properties of paracetamol are improved through the strategy of cocrystal formation.
    Graphical Abstract taken from @karkiImprovingMechanicalProperties2009.
  ],
)

#box(
  stroke: 2pt + green,
  inset: 1mm,
  [
    While many studies have focused on the role of dispersion, 25, 29, 30 few have critically assessed the shortcomings of semi-local functionals in detail.
    For example, the de-localisation or selfinteraction errors in DFT31 can often have a signiﬁcant effect on hydrogen-bonded systems. @reillyUnderstandingRoleVibrations2013
  ],
)

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

The harmonic treatment ignores effects due to anharmonic thermal motion and cell expansion. @reillyUnderstandingRoleVibrations2013

In the harmonic approximation, supercell phonon calculations show signiﬁcant deviations from the widely used Dulong-Petit law, as noted elsewhere. @reillyUnderstandingRoleVibrations2013

Assessing predictions of lattice energies requires careful consideration of vibrational, many-body dispersion and exact-exchange contributions. @reillyUnderstandingRoleVibrations2013

Exact exchange, which is rarely considered in @dft studies of molecular crystals, is shown to have a signiﬁcant contribution to lattice energies, systematically improving agreement between theory and experiment. @reillyUnderstandingRoleVibrations2013

Hybrid functionals are often not used in the study of cohesive properties of molecular crystals,
@al-saidiAssessmentVdWTSMethod2012, @otero-de-la-rozaBenchmarkNoncovalentInteractions2012
largely due to their additional computational cost, which particularly in a plane-wave basis can reach more than an order of magnitude larger than the corresponding semi-local functional. @reillyUnderstandingRoleVibrations2013

#box(
  stroke: 2pt + red,
  inset: 1mm,
  [The pharmaceutical industry spends considerable resources on high-throughput crystallization experiments to screen for polymorphs, // TODO citazione a S. L. Morissette et al., High-throughput crystallization: Polymorphs, salts, co-crystals and solvates of pharmaceutical solids. Adv. Drug Deliv. Rev. 56, 275–300 (2004).
    into which the target structure may decay.
    However, crystallization experiments do not probe thermodynamic stability, and conclusive studies of the impact of temperature changes after crystallization on the stability of polymorphs (i.e. their monotropic or enantiotropic nature) // TODO cite E. H. Lee, A practical guide to pharmaceutical polymorph screening & selection. Asian J. Pharm. Sci. 9, 163–175 (2014)
    are often prevented by limited sample quantities.
    Hence, there is the appeal of theoretical @csp // TODO cite 14 S. L. Price, Predicting crystal structures of organic compounds. Chem. Soc. Rev. 43,2098–2111 (2014).
    based on the thermodynamic stability, which promises to complement crystallization experiments // TODO cite 15 J. Nyman, S. M. Reutzel-Edens, Crystal structure prediction is changing from basic science to applied technology. Faraday Discuss. 211, 459–476 (2018).
    by exhaustively searching for competing polymorphs.
    @kapilCompleteDescriptionThermodynamic2022
  ],
)

#box(
  stroke: blue + 2pt,
  inset: 1mm,
  [
    The size and sign of nuclear quantum effects, anharmonicity, and cell expansion and flexibility, depend entirely on the compound and the polymorphs at hand, highlighting that rigorous @qti is indispensable for predicting phase stabilities and that molecular crystals are typically stabilized by a nontrivial interplay of different physical effects, whose individual importance is belied by the subtle resultant free energy differences. @kapilCompleteDescriptionThermodynamic2022
  ],
)

=== Machine Learning potentials

The facile generation of machine-learning potentials for a diverse set of polymorphic compounds—benzene, glycine, and succinic acid—and predictions of thermodynamic stabilities in qualitative and quantitative agreement with experiments highlight that predictive thermodynamic studies of industrially relevant molecular materials are no longer a daunting task. @kapilCompleteDescriptionThermodynamic2022

#box(
  stroke: 2pt + red,
  inset: 1mm,
  [
    Modern #glspl("mlp") // TODO cite 31 V. L. Deringer, M. A. Caro, G. Csányi, Machine learning interatomic potentials as emerging tools for materials science. Adv. Mater. 31, e1902765 (2019).
    permit accurately reproducing ab initio #glspl("pes") and dramatically reduce the cost of performing simulations approaching ab initio accuracy. // TODO cite 32 J. Lan et al., Simulating the ghost: Quantum dynamics of the solvated electron. Nat. Commun. 12, 766 (2021). Nat. Commun.12
    @kapilCompleteDescriptionThermodynamic2022
  ],
)

@kapilCompleteDescriptionThermodynamic2022 puts the cost of calculations into perspective,
highlighting that mixed ab initio and @mlp calculations require between three and ten times less the core hours of purely physical functionals.

== Research question

=== Ice polymorphs

Previous works combine @pi approaches with #glspl("mlp") for ice polymorphs
// TODO kapil 29 : V. Kapil, E. Engel, M. Rossi, M. Ceriotti, Assessment of approximate methods for anharmonic free energies. J. Chem. Theory Comput. 15, 5845–5857 (2019).
@chengInitioThermodynamicsLiquid2019

#set quote(block: true)

Taking into account the disruptive innovation brought in the molecular dynamics field by machine learning force fields, this thesis will try to answer the following question:
#quote[what is the performance of machine learning force fields in the simulation of water?]
The following chapters will introduce the theoretical foundations and the tools needed to pursue this question.

= Theory

== Phonons
The _dynamical matrix_:
$
  D(arrow(q)) :=
  1 / m sum_j e^(i arrow(q) dot arrow(r)_j^0)
  Phi(arrow(w)_i^0 - arrow(r)_j^0) dot arrow(epsilon)
$ <eq:dynamical-matrix>
$
  omega^2 arrow(epsilon) = D(arrow(q)) dot epsilon
$
$
  F(V,T) = U_0(V) + k_B T sum_(s=1)^3 sum_(arrow(q)) ln (
  planck.reduce omega_(arrow(q),s)(V)) / ( k_B T
  )
$
=== The small displacement method

#box(
  stroke: 2pt + red,
  inset: 1pt,
  [
    Let us consider the crystal in the ground state and displace one particle from its equilibrium position.
    Let the displacement be along the $x$ axis, $arrow(u) := (u, 0, 0)$.
    The forces acting on all particles in the system, including the displaced one, are given by:
  ],
)

$
  (f_i^x, f_i^y, f_i^z)
  = -Phi(arrow(r)_i^0) dot (u, 0, 0)
  = -u vec(phi_(0i)^(x x), phi_(0i)^(y x), phi_(0i)^(z x))
$

#box(
  stroke: 2pt + red,
  inset: 1pt,
  [
    from which we obtain the first columnt of the force constant matrix:
  ],
)

$
  phi_(0i)^(x x) = - f_i^x / u, quad
  phi_(0i)^(y x) = - f_i^y / u, quad
  phi_(0i)^(z x) = - f_i^z / u
$

To obtain the other columns, we need to repeat the procedure with displacements $(0,u,0)$ and $(0,0,u)$.
A set of three forces calculations is sufficient to compute the full force constants matrix.

In practical calculations, it is not possible to include all the infinite particles of a lattice;
therefore, we use a procedure that approximates the evaluation of the force constants matrix, with an error that can be systematically reduced.
We consider a _supercell_, which is an integer multiple of the primitive cell in the three spatial directions. We can define three multiplicative integers, one for each direction, $i,j,k$.
Named $arrow(a)_1, arrow(a)_2, arrow(a)_3$ the lattice vectors of the primitive cell,
the periodicity of the supercell is described by the multiplied lattice vectors:
$
  arrow(A)_1 := i arrow(a)_1, quad
  arrow(A)_2 := j arrow(a)_2, quad
  arrow(A)_3 := k arrow(a)_3
$
The reciprocal lattice vectors of the supercell are defined as follow:
$
  arrow(B)_1 := 2 pi (arrow(A)_2 times arrow(A)_3) / V, quad
  arrow(B)_2 := 2 pi (arrow(A)_3 times arrow(A)_1) / V, quad
  arrow(B)_3 := 2 pi (arrow(A)_1 times arrow(A)_2) / V,
$
with $V := (arrow(A)_1 times arrow(A)_2) dot arrow(A)_3$ as the volume of the supercell.
Following the periodicity of the finite supercell,
the displacement by an amount $arrow(u)$ of one particle
causes the displacement of all its periodic images too,
located at positions
$arrow(L) := n arrow(A)_1 + m arrow(A)_2 + l arrow(A)_3$,
with $n,m,l$ any three integers.
The obtained forces are:
$
  arrow(f)_i
  = - sum_(arrow(L)) Phi (arrow(r)_i^0 + arrow(L)) dot arrow(u)
$
where the sum runs over all possible $arrow(L)$ vectors distancing different supercells.
We do not have direct access to the elements of $Phi(arrow(r)_i^0)$,
but only to the force constants matrix of the supercell, defined as:
$
  Phi_"SC" (arrow(r)_i^0)
  := sum_(arrow(L)) Phi (arrow(r)_i^0 + arrow(L)).
$ <eq:force-constants-matrix-supercell>
This also implies that it is impossible to obtain the exact dynamical matrix,
which requires knowledge of every element of the force constants matrix.
We can rewrite @eq:dynamical-matrix as:
$
  D(arrow(q))
  &= 1 / m sum_j sum_arrow(L) e^(i arrow(q) dot (
    arrow(r)_j^0 + arrow(L)
  )) Phi(arrow(r)_j^0 + arrow(L)) \
  &= 1 / m sum_j e^(i arrow(q) dot arrow(r)_j^0)
  sum_arrow(L) e^(i arrow(q) dot arrow(L))
  Phi(arrow(r)_j^0 + arrow(L))
$
with the sum over $j$ running on the lattice sites $arrow(r)_j^0$ belonging to the supercell.
Consider now the term $sum_arrow(L) e^(i arrow(q) dot arrow(L)) Phi(arrow(r)_j^0 + arrow(L))$ and compare it with @eq:force-constants-matrix-supercell.
For every $arrow(q)$ such that $e^(i arrow(q) dot arrow(L)) = 1$ we obtain:

// TODO definisci vettori del reticolo
// TODO definisci matrice dinamica
// TODO finisci teoria

$
  D(arrow(q))
  = 1 / m sum_j e^(i arrow(q) dot arrow(r)_j^0) Phi_("SC")(arrow(r)_j^0)
$

Since $Phi_"SC"$ can be calculated,
at these particular $arrow(q)$ vectors the dynamical matrix is exact.
We can carachterize these special $arrow(q)$ vectors as the linear combinations of the reciprocal vectors of the supercell with integer coefficients $n',m',l'$.
We can verify this statement calculating, given $arrow(q) := n' arrow(B)_1 + m' arrow(B)_2 + l' arrow(B)_3$
and $arrow(L) := n arrow(A)_1 + m arrow(A)_2 + l arrow(A)_3$,
$arrow(q) dot arrow(L) = (n n' + m m' + l l') 2pi$,
which yields $e^(i arrow(q) + arrow(L)) = 1$.

Increasing the size of the supercell, the Brillouin zone is populated with more exact phonons;
inbetween those points we obtain a _Fourier interpolation_,
which becomes more and more accurate as we increase the size of the supercell.
Eventually, the supercell is so large that the force constants matrix is negligible at its edges;
when this happens, the only term contributing in the sum in @eq:force-constants-matrix-supercell is that with $arrow(L) = arrow(0)$;
as we increase the size of the supercell, the supercell force constants matrix asymptotically approaches the force constants matrix,
$Phi_"SC" (arrow(r)_i^0) tilde.eq Phi(arrow(r)_i^0)$.
In this limit, the Fourier interpolation is accurate everywhere.

== DFT
Kohn-Sham equations:
$
  [-1 / 2 nabla^2 + v_"KS" [rho] (arrow(r))] psi_n (arrow(r))
  = epsilon_n psi_n (arrow(r))
$
The aim is to find the lowest $N/2$ eigenstates of the Hamiltonian.
These are of the type:
$
  psi_n = sum_(i=1)^L c_i^n phi_i
$
The energy can be obtained from:
$
  E = 2 sum_(n=1)^(N / 2) epsilon_n - 1 / 2 integral_V (rho (
    arrow(r)'
  ) rho(arrow(r))) / (|arrow(r) - arrow(r)'|) dif^3 arrow(r)' dif^3 arrow(r)
  - integral_V (delta E_"xc") / (delta rho(arrow(r))) rho(arrow(r)) dif^3 arrow(r) + E_"xc"
$
== Molecular dynamics
=== The Verlet algorithm
The Verlet algorithm is a technique to generate the trajectory of interacting particles obeying the Newton's equations of motion. @alfeNotesStatisticalComputational2023
It is a discretization of Newton's equations of motion:
$
  arrow(f)_i = M dot.double(arrow(r))_i = - (partial U(
    {arrow(r)}
  )) / (partial arrow(r)_i)
$ <eq:verlet-newton>
Let us consider the Taylor expansion of the position of particle $i$ at time $t$, $arrow(r)_i (t)$, computed with forward and backward differences:
$
  arrow(r)_i (t + delta t)
  = arrow(r)_i (t) + dot(arrow(r))_i (
    t
  ) delta t + 1 / 2 dot.double(arrow(r))_i (t) (
    delta t
  )^2 + 1 / (3!) dot.triple(arrow(r))_i (t) (delta t)^3 + O((delta t)^4)
$
$
  arrow(r)_i (t - delta t)
  = arrow(r)_i (t) - dot(arrow(r))_i (
    t
  ) delta t + 1 / 2 dot.double(arrow(r))_i (t) (
    delta t
  )^2 - 1 / (3!) dot.triple(arrow(r))_i (t) (delta t)^3 + O((delta t)^4)
$
Summing the two equations side by side, we obtain:
$
  arrow(r)_i (t + delta t) + arrow(r)_i (t - delta t)
  = 2 arrow(r)_i + dot.double(arrow(r))_i (t) (delta t)^2 + O((delta t)^4)
$
Consider the expression of $dot.double(arrow(r))_i$ in terms of $arrow(f)_i$ from @eq:verlet-newton, $dot.double(arrow(r))_i (t) = 1/M arrow(f)_i (t)$; substituting, we obtain:
$
  arrow(r)_i (t + delta t)
  = 2 arrow(r)_i (t) - arrow(r)_i (t - delta t) + 1 / M arrow(f)_i (t) (
    delta t
  )^2 + O((delta t)^4)
$ <eq:verlet-algorithm>
@eq:verlet-algorithm is known ad the Verlet algorithm.
@verletComputerExperimentsClassical1967[Eq. (4)]

We can re-express the equation in terms of the velocities
$
  arrow(v)_i (t) = (arrow(r)_i (t + delta t) - arrow(r)_i (
    t - delta t
  )) / (2 delta t)
$
$
  - arrow(r)_i (t - delta t)
  = arrow(v)_i (t) dot 2 delta t
  - arrow(r)_i (t + delta t)
$
$
  2 arrow(r)_i (t + delta t)
  = 2 arrow(r)_i (t)
  + 2 arrow(v)_i (t) delta t
  + 1 / M arrow(f)_i (t) (delta t)^2
  + O((delta t)^4)
$
$
  arrow(r)_i (t + delta t)
  = arrow(r)_i (t)
  + arrow(v)_i (t) delta t
  + 1 / (2M) arrow(f)_i (t) (delta t)^2
  + O((delta t)^4)
$
This expression gives access to the positions at time $t + delta t$ with just the knowledge of positions, velocities, and forces at time $t$.

== Neural networks

= Results I: model assessment

In this chapter we test the MACE calculators on small, known systems,
about which literature is abundant and a direct comparison of results is possible.
We analyze the geometrical and vibrational features of the water molecule and the water dimer after optimization with each calculator.
The binding energy of the dimer is calculated and compared with reference @dft methods.
This setting already allows to quantify the benefits of a fine-tuned @mlp like MACE-ICE13-1 compared with MACE foundation models.

// Seconda parte: simulazione per qualche sistema standard in cui l'approccio analitico funziona bene

// Terza parte: simulazione per qualche sistema in cui l'approccio analitico non funziona bene

== The water molecule <sec-molecule>

The first task is the optimization of the geometry and calculation of vibrational properties of the water molecule.
We studied the relaxation of the geometry of a single water molecule, also referred to as the monomer, and compared the results with physical values of the gas phase molecule. The optimization is tackled with two concurring methods:
static local *minimization of the potential energy* and *analysis of vibrational properties* to assess the dynamical stability.

#figure(
  image("simulazioni/02_water/01_molecule/MACE-ICE13-1/final.png", width: 70%),
  caption: [
    Render of the geometry of the water molecule,
    optimized using MACE-ICE13-1.
  ],
)

=== Geometry optimization

To rapidly put to the test the calculators, geometrical values of the relaxed configuration have been computed and compared with references in literature.

The first value concerns the charachteristic $#ce("HOH")$ *bend angle* of the molecule;
the accurate description of this physical value is a required test to ensure the validity of the models.
As can be seen in @table:hoh-angle, the MACE-MP-0 large model most accurately describes the bend angle of the molecule, while the small model is the most distant from reference @eisenbergStructurePropertiesWater2005 @PhysicalChemistryWater2020.

The second value is the $#ce("OH")$ *bond length*.
@table:oh-bond-length shows that MACE-ICE13-1 gives the most accurate description according to reference @eisenbergStructurePropertiesWater2005 @PhysicalChemistryWater2020, and MACE-MP-0 small again is less accurate than the other models.

#let monomer_angles_table = csv("simulazioni/02_water/01_molecule/angles.csv")
#figure(
  table(
    columns: monomer_angles_table.first().len(),
    table.header([Model], [$#ce("HOH")$ Angle (°)], [Discrepancy (°)]),
    ..monomer_angles_table.slice(1).flatten(),
  ),
  caption: [
    The calculated $#ce("HOH")$ bend angle for the gas phase water monomer, and difference with respect to reference.
    small, medium and large refer to MACE-MP-0.
  ],
) <table:hoh-angle>

#figure(
  image(
    "simulazioni/02_water/01_molecule/Grafici/angle_convergence_mace_mp_0_large.svg",
    width: 90%,
  ),
  caption: [Convergence of the $#ce("HOH")$ angle with respect to the `fmax` parameter.],
)

#let monomer_bond_length_table = csv("simulazioni/02_water/01_molecule/bond_lengths.csv")
#figure(
  table(
    columns: monomer_bond_length_table.first().len(),
    table.header(
      [Model],
      [Bond length ($angstrom$)],
      [Discrepancy ($angstrom$)],
    ),
    ..monomer_bond_length_table.slice(1).flatten(),
  ),
  caption: [
    Calculated $#ce("OH")$ bond lengths for the gas phase water monomer, and difference with respect to reference.
  ],
) <table:oh-bond-length>

The following sections detail the procedure to assess the convergence and stability of the obtained geometries, as well as explore the energetics of the molecule.

=== The optimization algorithm

The procedure to optimize the geometry of the molecule requires the correct setting of a few parameters, namely:

- the calculator and the data type
- the starting geometry and cell properties
- the optimizer algorithm
- the force threshold to stop the optimization
- the maximum number of steps of the optimization

The calculators chosen for the optimization are MACE-MP-0, in the small, medium and large variants, and MACE-ICE13-1.
As advised by the calculators, the float64 data type was chosen for optimization, where float32 is advised for molecular dynamics runs for improved speed.

The initial geometry was set up so that the HOH atoms formed a right angle, and the OH distance was $1 angstrom$.
By dealing with a force field, it was possible to set up a cell in vacuum without periodic boundary conditions.

@ase offers a variety of different optimisers, from which we have chosen BFGS, a local optimisation algorithm of the quasi-Newton category, where the forces of consecutive steps are used to dynamically update a Hessian describing the curvature of the potential energy landscape.
The optimiser accepts two important input parameters. The first if _fmax_, the force threshold, defined in $"eV" angstrom^(-1)$.
The convergence criterion is that the force on all individual atoms should be less than fmax:
#footnote[https://wiki.fysik.dtu.dk/ase/ase/optimize.html]
$
  max_a |arrow(F)_a| < f_"max"
$

For the present purposes, `fmax=1e-8` yielded satisfactory results for the different calculator models tested.

The second parameter is the number of steps after which to stop the optimization procedure.
If the steps employed to optimize the geometry, given the desired fmax, are larger than the steps parameter, the procedure is halted and the geometry is considered as not converged.
A value of `steps=1000` is largely above the actual number of steps required for convergence of the geometry.

#[
  #show figure: set align(left)
  #figure(
    ```python
    from ase.optimize import BFGS
    opt = BFGS(
      atoms,
      logfile="optimization.log",
      trajectory="optimization.traj"
    )
    opt.run(fmax=1e-8, steps=1000)
    ```,
    caption: [
      Initialisation and run of the optimizer.
    ],
  )
]

=== Molecular vibrations and assessment of the dynamical stability

The nuclei of molecules, far from occupying fixed positions with respect to each other, are in continual state of vibration, even at 0°K.
An important feature of these vibrations is that they can be described by a limited number of basic vibrations known as the normal modes.
A normal mode is a vibration in which all the nuclei oscillate with the same frequency and the same phase.
The water molecule has three normal modes and every possible vibration of the molecule can be described as a superposition of these three modes.

The normal modes of vibration of water are shown in @fig:h2o-normal-modes.
Because the motion of the nuclei in the $nu_1$ and $nu_3$ vibrations is nearly along the direction of the $#ce("O - H")$, these modes are often referred to as $#ce("O - H")$ stretching vibrations.
Similarly, because the $#ce("H")$ nuclei in $nu_2$ move in directions almost perpendicular to the bonds, $nu_2$ is referred to as the $#ce("H - O - H")$ bending vibration.
In fact, $nu_1$ involves a small amount of $#ce("H - O - H")$ bending, and $nu_2$ involves a small amount of $#ce("O - H")$ stretching.
The mode $nu_3$ is called the asymmetric stretching vibration to distinguish it from the symmetric stretching vibration $nu_1$. @eisenbergStructurePropertiesWater2005[§1.1 (d)]

#figure(
  image("thesis/imgs/eisenberg2005_fig1.1.gif", width: 50%),
  caption: [
    The normal modes of vibration of $#ce("H2O")$.
    Taken from @eisenbergStructurePropertiesWater2005[Fig. 1.1].
  ],
) <fig:h2o-normal-modes>

In this work, the vibration modes were studied under the harmonic approximation.
Let us denote by $G(v_1, v_2, v_3)$ the energy above the vibrationless equilibrium state of the state with quantum numbers $v_1, v_2$ and $v_3$. Then
$
  G(v_1, v_2, v_3)
  = sum_(i=1)^3 omega_i (v_i + 1 / 2)
  + sum_(i=1)^3 sum_(k >= i)^3 x_(i k) (v_i + 1 / 2) (v_k + 1 / 2)
$
where the sums are over normal modes.
The $omega$s in this equation are often called the _harmonic frequencies_;
they are the frequencies with which the molecule would vibrate if its vibrations were perfectly harmonic.
The $x$s are the _anharmonic constants_ and describe the effect on the vibrational frequencies of the departure from purely harmonic form of the vibrations.
@table:molecule-omega and @table:molecule-omega-errors
contain the harmonic frequencies for $#ce("H2O")$ and the errors for each model compared with reference data.

#large_box(
  grid(
    columns: 2,
    gutter: 10pt,
    [
      #let molecule_omega_table = csv("simulazioni/02_water/01_molecule/Analysis/omega.csv")
      #figure(
        table(
          columns: molecule_omega_table.first().len(),
          table.header(
      // ..molecule_omega_table.first(),
      [Model], $omega_1$, $omega_2$, $omega_3$
    ),
          ..molecule_omega_table.slice(1).flatten()
        ),
        caption: [
          The harmonic frequencies of the normal modes of the water molecule,
          expressed in $"cm"^(-1)$.
          Reference data from @eisenbergStructurePropertiesWater2005[Table 1.4].
        ],
      ) <table:molecule-omega>
    ],
    [
      #let molecule_omega_errors_table = csv("simulazioni/02_water/01_molecule/Analysis/omega_errors.csv")
      #figure(
        table(
          columns: molecule_omega_errors_table.first().len(),
          table.header(
      // ..molecule_omega_table.first(),
      [Model], $omega_1$, $omega_2$, $omega_3$
    ),
          ..molecule_omega_errors_table.slice(1).flatten()
        ),
        caption: [
          Errors on the harmonic frequencies of the normal modes of the water molecule,
          expressed in $"cm"^(-1)$.
        ],
      ) <table:molecule-omega-errors>
    ],
  ),
)

Studying the vibrational properties of the geometry obtained at the end of the optimization procedure also allows us to assess if the final geometry is a stable or unstable configuration.
The vibrational modes are calculated from a finite difference approximation of the Hessian matrix, displacing atoms according to a parameter named `delta`, measured in $angstrom$.
#footnote[https://wiki.fysik.dtu.dk/ase/ase/vibrations/modes.html]

The following figures represent the frequencies of the vibration modes of the optimized water molecule obtained with @ase plotted against the displacement parameter `delta`, using different calculator models and values of the `fmax` parameter.
Frequencies are indexed in ascending order, and *imaginary frequencies*, representing unstable configurations, are shown as negative values in the graphs.
Inclusion of dispersion contributions in calculators leads to minimal differences in the
converged configurations.

The stability of configurations is dependent in a discriminant way on the `delta` and `fmax` parameters.
The analysis confirms that the value of `fmax=1e-4` yields unstable configurations, while `fmax=1e-8` is appropriate to obtain stable configurations.
Moreover, the displacement `delta` shall be smaller than `1e-4` to ensure convergence of calculations.

#large_figure(
  grid(
    columns: 3,
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small fmax=1e-1.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small fmax=1e-4.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small fmax=1e-8.svg"),
  ),
  caption: "Convergence of the small model with respect to fmax",
)

#large_figure(
  grid(
    columns: 3,
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium fmax=1e-1.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium fmax=1e-4.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium fmax=1e-8.svg"),
  ),
  caption: "Convergence of the medium model with respect to fmax",
)

#large_figure(
  grid(
    columns: 3,
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large fmax=1e-1.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large fmax=1e-4.svg"),
    image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large fmax=1e-8.svg"),
  ),
  caption: "Convergence of the large model with respect to fmax",
)

// #large_figure(
//   grid(
//     columns: 3,
//     image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small D fmax=1e-1.svg"),
//     image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small D fmax=1e-4.svg"),
//     image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 small D fmax=1e-8.svg"),
//   ),
//   caption: "Convergence of the small model with dispersion with respect to fmax",
// )

// #large_figure(
//   grid(
//     columns: 3,
//     image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium D fmax=1e-1.svg"),
//     image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium D fmax=1e-4.svg"),
//     image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium D fmax=1e-8.svg"),
//   ),
//   caption: "Convergence of the medium model with dispersion with respect to fmax",
// )

// #large_figure(
//   grid(
//     columns: 3,
//     image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large D fmax=1e-1.svg"),
//     image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large D fmax=1e-4.svg"),
//     image("simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 large D fmax=1e-8.svg"),
//   ),
//   caption: "Convergence of the large model with dispersion with respect to fmax",
// )

Adding the dispersion correction to the models in this step produces insignificant differences in the results, so their graphs are omitted for brevity.

==== MACE-ICE13-1

The MACE-ICE13-1 model is fine-tuned from the MACE-MP-0 medium model with
dispersion.

Convergence of results for the MACE-ICE13-1 model is achieved for displacements below $10^(-4) angstrom$.
The imaginary frequency observable in @fig-monomer-vibrations-mace-ice13-1 corresponds to negligible energy, as can be observed in the representative output in @code-monomer-vibrations-mace-ice13-1-output and therefore does not pose a issue.

#large_box(
  grid(
    columns: 2,
    gutter: 7pt,
    [#figure(
        image("simulazioni/02_water/01_molecule/Grafici/MACE-ICE13-1 fmax=1e-8.svg"),
        caption: [Frequencies calculated with the MACE-ICE13-1 model with the most restrictive `fmax`.],
      ) <fig-monomer-vibrations-mace-ice13-1>],
    [
      #figure(
        raw(
          read("simulazioni/02_water/01_molecule/MACE-ICE13-1/H2O_delta=1e-06_summary.txt"),
        ),
        caption: [
          Vibrational analysis representative output from ASE for the normal modes of the water molecule, computed using MACE-ICE13-1.
        ],
      ) <code-monomer-vibrations-mace-ice13-1-output>
    ],
  ),
)


#figure(
  ```python
  for d in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]:
      vib = Vibrations(
        atoms,
        delta=d,
        nfree=4,
        name=f"vib_delta={d}"
      )
      vib.run()
      vib.summary(log=f"H2O_delta={d}_summary.txt")
      vib.clean()
  ```,
  caption: [
    Computation of vibrational properties for different values of the displacement of atoms.
  ],
)

==== Zero-point vibrational energy

The @zpe from the computation models is shown in @table:zpe.
The MACE-ICE13-1 model shows the best agreement with reference data from literature @eisenbergStructurePropertiesWater2005, with a discrepancy of $0.01 "eV"$.

#let zero_point_energies_table = csv("simulazioni/02_water/01_molecule/zero_point_energies.csv")
#figure(
  table(
    columns: zero_point_energies_table.first().len(),
    // table.header(..zero_point_energies_table.first()),
    table.header([Model], [@zpe ($"eV"$)], [Discrepancy (eV)]),
    ..zero_point_energies_table.slice(1).flatten(),
  ),
  caption: [
    Zero-point energies for the employed calculators and reference value @eisenbergStructurePropertiesWater2005.
    The difference between calculated and reference value is also shown in the last column.
    small, medium and large models refer to MACE-MP-0.
  ],
) <table:zpe>

The results so far indicate that the medium, large MACE-MP-0, and MACE-ICE13-1 approximate better the properties of the water molecule,
while the small MACE-MP-0 model is the worst of them in this regard.

== The water dimer

After analyzing the properties of the water molecule it is time to study a slightly more elaborate system.
The water dimer exposes a bigger number of geometric features and vibrational modes.
Moreover, a charactheristic physical quantity, the binding energy, allows to quantitatively evaluate the performance of the models.
Resources on the water dimer are abundant, on both the theory and experiment sides.
@mukhopadhyayWaterDimerII2018
@barnettBornOppenheimerMoleculardynamicsSimulations1993
@mukhopadhyayWaterDimerExperimental2015
@schutzWaterDimerInteraction1997
@klopperComputationalDeterminationEquilibrium2000
@hobzaReliableTheoreticalTreatment1999
Harmonic approximation reference values @kalesckyLocalVibrationalModes2012 were used for comparisons in @sec:dimer-vibrations.

=== Geometry optimization

The optimization was performed starting from an initial configuration
approximately reproducing the geometry represented in @dimer-structure. The BFGS
optimizer was used with the three MACE-MP-0 models and the MACE-ICE13-1 model, with a
force threshold of `fmax=1e-8` (units $"eV"slash angstrom$).
The code describing the initial geometry is available in @init.xyz.
The results comparing values and errors of the soft, intermolecular, features of the optimized geometry with respect to literature
reference are available in @dimer-geometry-table and @dimer-geometry-errors.

#large_box(
  grid(
    columns: (1fr, 1.5fr),
    column-gutter: 10pt,
    [
      #figure(
        image("thesis/imgs/klopper-fig1.gif"),
        caption: [
          The equilibrium structure of the water dimer. (Image taken from
          @klopperComputationalDeterminationEquilibrium2000)
        ],
      ) <dimer-structure>
    ],
    [
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
          The initial geometry for the optimization of the water dimer,
          `init.xyz`.
        ],
      ) <init.xyz>
    ],
  ),
)

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

The overall results exhibit a better accuracy of the MACE-MP-0 large and MACE-ICE13-1 models,
holding first and second place with lowest discrepancies with respect to reference most of the time,
while the small and medium MACE-MP-0 models lag behind.

A discussion on the importance of parameters can be made,
considering that $alpha$ is an H-bonded geometrical parameter,
and can be deemed more important in the evaluation of the performance of the models,
with respect to $beta$ that is not H-bonded,
and is allowed a relatively bigger range of motion.

#figure(
  image("simulazioni/02_water/02_dimer/01_optimize/MACE-ICE13-1/final.png"),
  caption: [
    Render of the final geometry of the water dimer,
    optimized using MACE-ICE13-1.
  ],
)

=== Vibrations analysis <sec:dimer-vibrations>

Similarly to the procedure adopted for the monomer,
normal modes of the water dimer were analyzed to assess the dynamical stability of the optimized geometry.
A displacement smaller than $10^(-4) angstrom$ is required for all the models to converge to stable configurations.

#large_figure(
  grid(
    columns: 2,
    image("simulazioni/02_water/02_dimer/01_optimize/Grafici/small.png"),
    image("simulazioni/02_water/02_dimer/01_optimize/Grafici/medium.png"),

    image("simulazioni/02_water/02_dimer/01_optimize/Grafici/large.png"),
    image("simulazioni/02_water/02_dimer/01_optimize/Grafici/MACE-ICE13-1.png"),
  ),
  caption: [
    Structure optimization of the water dimer and vibration analysis.
    Negative values represent imaginary frequencies.
  ],
)

@table:dimer-harmonic-frequencies-comparison and @table:dimer-harmonic-frequencies-errors show the comparison and the errors of the harmonic frequencies of the 12 normal modes of the water dimer with respect to reference @kalesckyLocalVibrationalModes2012;
reference values were calculated using @ccsdt with a @cbs.
Deviations from reference are also graphed in @fig:dimer-harmonic-frequencies-errors.

#let dimer_harmonic_frequencies_comparison = csv("simulazioni/02_water/02_dimer/01_optimize/Analisi/harmonic_frequencies_comparison.csv")
#figure(
  table(
    columns: dimer_harmonic_frequencies_comparison.first().len(),
    table.header(..dimer_harmonic_frequencies_comparison.first()),
    ..dimer_harmonic_frequencies_comparison.slice(1).flatten()
  ),
  caption: [Comparison of the harmonic normal modes frequencies obtained through MACE calculators and reference @kalesckyLocalVibrationalModes2012. Frequency units are in $"cm"^(-1)$.],
) <table:dimer-harmonic-frequencies-comparison>

#let dimer_harmonic_frequencies_errors = csv("simulazioni/02_water/02_dimer/01_optimize/Analisi/harmonic_frequencies_errors.csv")
#figure(
  table(
    columns: dimer_harmonic_frequencies_errors.first().len(),
    table.header(..dimer_harmonic_frequencies_errors.first()),
    ..dimer_harmonic_frequencies_errors.slice(1).flatten()
  ),
  caption: [Difference between the harmonic normal modes frequencies obtained through MACE calculators and reference @kalesckyLocalVibrationalModes2012. Frequency units are in $"cm"^(-1)$.],
) <table:dimer-harmonic-frequencies-errors>

#figure(
  image("/simulazioni/02_water/02_dimer/01_optimize/Grafici/harmonic_frequencies_errors_barchart.svg"),
  caption: [
    Deviation of the harmonic frequencies with respect to reference @kalesckyLocalVibrationalModes2012 for each MACE calculator.
  ],
) <fig:dimer-harmonic-frequencies-errors>

In @table:dimer-frequencies-sum-absolute-errors we make a summary of the predictive power of each model of the frequencies of the normal modes of the dimer, as compared to reference @kalesckyLocalVibrationalModes2012.
MACE-ICE13-1 demonstrates the best overall adherence to the harmonic frequencies.

#let dimer_harmonic_frequencies_mae = csv("simulazioni/02_water/02_dimer/01_optimize/Analisi/mae.csv")
#figure(
  table(
    columns: dimer_harmonic_frequencies_mae.first().len(),
    table.header([Model], [MAE ($"cm"^(-1)$)]),
    ..dimer_harmonic_frequencies_mae.slice(1).flatten()
  ),
  caption: [
    @mae computed as the difference between the value of the frequencies obtained with caluculator models and reference @kalesckyLocalVibrationalModes2012.
  ],
) <table:dimer-frequencies-sum-absolute-errors>

=== ZPE
#let dimer_zpe = csv("simulazioni/02_water/02_dimer/01_optimize/Analisi/zpe.csv")
#figure(
  table(
    columns: dimer_zpe.first().len(),
    table.header(..dimer_zpe.first()),
    ..dimer_zpe.slice(1).flatten()
  ),
  caption: [
    #glspl("zpe") of the water dimer in the harmonic approximation and comparison with reference @kalesckyLocalVibrationalModes2012. Energies are in units of eV.
  ],
)

=== Binding energy

The binding energy is calculated as

$ Delta E_2 := E_2 - 2 E_1 $

where $E_2$ is the energy of the dimer, $E_1$ is the energy of one molecule.

To determine the equilibrium geometry of the dimer from the study of the binding
energy, $Delta E_2$ must be minimized with respect to the internal coordinates.
To reduce the complexity of the task, the optimization is constrained to the
single degree of freedom of the $#ce("O - O")$ distance.
#footnote[
  The FixBondLength class from @ase was employed to enforce the constraint.
  See https://wiki.fysik.dtu.dk/ase/ase/constraints.html#the-fixbondlength-class
]
@klopperComputationalDeterminationEquilibrium2000
Distances between 2 and 6 $angstrom$ were sampled to build the graph,
with a denser sampling near the respective minima for each of the calculators.
Once the two molecules are spaced apart of the fixed distance, a relaxation is performed and the potential energy of the system is calculated.

In @fig:binding-energy the binding energy is graphed against the $#ce("O - O")$ distance, $r_"OO"$,
for MACE-MP-0 medium and MACE-ICE13-1.
The equilibrium distance according to @dykeStructureWaterDimer1977 @mukhopadhyayWaterDimerII2018, $r_"OO" = 2.98 angstrom$, is most accurately predicted by MACE-ICE13-1, while MACE-MP-0 medium is off by about a quarter of a $angstrom$.

Typical values for the binding energy obtained from @dft simulations at
equilibrium position range between -20 and -12 kJ/mol (approximately -0.2 to
-0.12 eV) @sprikInitioMolecularDynamics1996 @mukhopadhyayWaterDimerII2018.
Reference @curtissStudiesMolecularAssociation1979 reports the experimental binding energy for the water dimer to be $-5.44 plus.minus 0.7 "kcal/mol"$ that corresponds to $0.24 plus.minus 0.03 "eV/particle"$.
As MACE-ICE13-1 is trained on revPBE-D3, its predicted binding energy lies in the expected region according to @dft simulations.
See @fig:binding-energy-mukhopadhyay2018-fig5 for a representation of typical values of the binding energy obtained through @dft methods.

#large_box(
  grid(
    columns: 2,
    gutter: 10pt,
    [#figure(
        image("simulazioni/02_water/02_dimer/02_binding_energy/binding_energy.svg"),
        caption: [
          Dimer binding energy for different calculators. The vertical dashed line
          corresponds to the reference value for the equilibrium OO distance,
          corresponding to the experimental value of $2.98 angstrom$
          @dykeStructureWaterDimer1977 @mukhopadhyayWaterDimerII2018. The horizontal
          dashed line corresponds to $-0.24 "eV"$
          @curtissStudiesMolecularAssociation1979.
        ],
      ) <fig:binding-energy>
    ],
    [#figure(
        image("thesis/imgs/mukhopadhyayWaterDimerII2018_fig5.png"),
        caption: [
          Image from reference @mukhopadhyayWaterDimerII2018[Fig. 5],
          _The interaction energy of the water dimer for the equilibrium geometry configuration calculated at different center of mass separations between the two monomers_, using @dft and other lower level methods calculations.
        ],
      ) <fig:binding-energy-mukhopadhyay2018-fig5>
    ],
  ),
)

= Results II: crystal structures
#box(
  stroke: red + 2pt,
  inset: 1mm,
  [
    This chapter presents the results of the study, focusing on the properties of molecular crystals modeled using machine learning potentials.
    The analyses include geometry optimization, vibrational analysis, lattice energy computations, and crystal phonons calculations.
  ],
)

== Lattice energy

The work in this section is based on @dellapiaDMCICE13AmbientHigh2022.
The correct estimate of the lattice energies, absolute and relative, is of great interest for water and molecular crystals in general.
According to the article, the best values for the energies are consiered to be the @dmc calculations.
Regarding @dft methods, @guptaPhononsAnomalousThermal2018 found that the revPBE functional gives the closest structure and dynamical properties of ice Ih.

#figure(
  image("strutture/ICE13/Ih/Ih.png"),
  caption: [
    Render of the ice Ih cell.
  ],
)

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

$ E_"crystal" := E / N_(#ce("H2O")) $

The quantity $E_"gas"$ is calculated in the same manner as in @sec-molecule,
however, with the distinction that optimization was not performed;
to align the computed results with the reference paper, the same fixed gas phase molecule was used,
referenced as Patridge 1997 from the original authors,
to correctly reproduce the computational setup.

#figure(
  image("simulazioni/02_water/03_ICE13_lattice_energies/absolute_lattice_energy.svg"),
  caption: [
    Performance of MACE for the 13 ice polymorphs considered on the absolute lattice energy, compared with reference models from @dellapiaDMCICE13AmbientHigh2022.
  ],
)

The MACE-ICE13-1 model achieves a @mae with respect to DMC of $0.90 "kJ/mol"$.
@xc functionals are generally classified as good if their @mae is ≲4 kJ/mol (≲2 kJ/mol) for the absolute (relative) lattice energy. @dellapiaDMCICE13AmbientHigh2022[§III]

A possible improvement in the quality of calculations could be investigated employing bigger supercells for the calculation of lattice energies.
We could expect a different yield for the #glspl("mlp") as force field potentials are particularly suitable for dealing with a large number of particles, and the Ewald summation might be relevant for the quality of results in the present setting.

#large_figure(
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
@fig:relative-lattice-energy shows the computed relative lattice energy for the 13 selected polymorphs of water.
MACE-ICE13-1 shows a great adherence to its reference potential, revPBE-D3, positioning itself closely to the @dmc quality of results.

#figure(
  image("simulazioni/02_water/03_ICE13_lattice_energies/relative_lattice_energy.svg"),
  caption: [
    Performance of MACE for the 13 ice polymorphs considered on the relative lattice energy, compared with reference models. @dellapiaDMCICE13AmbientHigh2022
  ],
) <fig:relative-lattice-energy>

== Crystal phonons

The current implementation of phonons calculation in @ase
#footnote[https://wiki.fysik.dtu.dk/ase/ase/phonons.html]
is outclassed by Phonopy. #footnote[See https://gitlab.com/ase/ase-workshop-discussion/-/issues/7#note_245747917 and https://gitlab.com/ase/ase/-/issues/1235]

Normal modes of vibration are calculated using the so-called *small displacement method*.
This method is increasingly more accurate with bigger and bigger supercells.

The first thing to do is to build a supercell.
The Ih ice structure was analyzed;
it is composed of 36 atoms, that is 12 water molecules.
// TODO Definisci struttura del ghiaccio Ih
The biggest supercell that the GPU cuda version of MACE can handle is the $3 times 3 times 3$.
Using MACE on CPU allows us to employ also $4 times 4 times 4$ supercells, at the expense of a significantly increased computation time.

=== Band structure

To build the band structure of ice Ih,
we have to define the band path on which the frequencies are calculated.
The crystal form of ice Ih is hexagonal.
#footnote[https://en.wikipedia.org/wiki/Phases_of_ice#Known_phases]
The Brillouin zone of the hexagonal lattice is also hexagonal.
In the Brillouin zone one can find several points of high symmetry that are of special interest, called _critical points_.
Namely: #footnote[https://en.wikipedia.org/wiki/Brillouin_zone#Critical_points]

#grid(
  columns: 2,
  align: horizon,
  [
    - $Gamma$: center of the Brillouin zone
    - $A$: center of a hexagonal face
    - $H$: corner point
    - $K$: middle of an edge joining two rectangular faces
    - $L$: middle of an edge joining a hexagonal and a rectangular face
    - $M$: center of a rectangular face
  ],
  [
    #figure(
      image("strutture/ICE13/Ih/brillouin_zone_special_points.svg"),
      caption: [Special points of the Brillouin zone of ice Ih.],
    )
  ],
)

The bandpath is chosen following the reference article @guptaPhononsAnomalousThermal2018, that is sampling the Brillouin zone passing through critical points $Gamma, A, K, H, M, L, Gamma$ by means of straight line segments.
See @fig:bandpath-gupta for the graphical representation of the bandpath.
The line segments are sampled with a number of points, named q-points.
The number of q-points in each path including end points is chosen as 101.

#large_box(
  grid(
    columns: 2,
    align: horizon,
    gutter: 10pt,
    [#figure(
        image("simulazioni/02_water/04_crystal_phonons/phonopy/Grafici/bandstructure_mace-ice13-1_s3_gupta_full.svg"),
        caption: [
          The bandstructure of ice Ih, computed using MACE-ICE13-1 and a $3 times 3 times 3$ supercell.
        ],
      )],
    [#figure(
        image("simulazioni/02_water/04_crystal_phonons/phonopy/Grafici/bandpath_gupta.svg"),
        caption: [The bandpath chosen by @guptaPhononsAnomalousThermal2018.],
      ) <fig:bandpath-gupta>],

    [#figure(
        image("simulazioni/02_water/04_crystal_phonons/phonopy/Grafici/bandstructure_mace-ice13-1_s3_gupta.svg"),
        caption: [
          Phonon bandstructure of ice Ih, computed using MACE-ICE13-1.
        ],
      ) <fig:phonons-bandstructure-ice-ih-mace-ice13-1-zoom>],
    [#figure(
        image("thesis/imgs/gupta2018_fig4_H2O.png"),
        caption: [Phonon bandpath dispersion reference, taken from @guptaPhononsAnomalousThermal2018[Fig. 4].],
      ) <fig:phonons-bandstructure-ice-ih-gupta>],
  ),
)

The bandstructure calculations result in a approximate reproduction of reference data, as can be observed comparing @fig:phonons-bandstructure-ice-ih-mace-ice13-1-zoom and @fig:phonons-bandstructure-ice-ih-gupta.
Frequencies are significantly higher than reference; the source of this discrepancy is not clear, and might be due to improper volume optimization of our geometry or a built-in issue of the calculator; indeed, a compression of the material results in higher overall frequencies.
This issue has yet to be investigated at the time of writing.
Frequencies calculated with MACE-MP-0, shown in @fig:phonons-bandstructure-ice-ih-mace-mp-0-zoom, exhibit even higher frequencies and are reputed as lower quality for the current analysis.
A visual qualitative analysis of the band structure produced with MACE-MP-0 exposes several different behaviours, particularly at zone boundary points;
most notably, see the disalignment of bands around point $M$, and at point $A$.

As reference @guptaPhononsAnomalousThermal2018 employed a $2 times 2 times 2$ supercell for its calculations, we tested the convergence of calculations with respect to a higher supercell.
Limitation of resources allowed the calculation of the frequencies using at most $3 times 3 times 3$ supercell.
The comparison of the results with the different supercells is shown in @fig:phonons-bandstructure-ice-ih-mace-ice13-1-compare-s2-s3.

#large_box(
  grid(
    columns: 2,
    gutter: 10pt,
    [#figure(
        image("simulazioni/02_water/04_crystal_phonons/phonopy/Grafici/bandstructure_mace-ice13-1_s2-s3_gupta_zoom.svg"),
        caption: [Comparison of the bandstructures computed with $2 times 2 times 2$ and $3 times 3 times 3$ supercells.],
      ) <fig:phonons-bandstructure-ice-ih-mace-ice13-1-compare-s2-s3>],
    [
      #let original = read("simulazioni/02_water/04_crystal_phonons/phonopy/Grafici/bandstructure_mace-mp-0_s2_gupta.svg")
      #let changed = original.replace(
        "#ff0000",
        blue.to-hex(),
      )
      #figure(
        image.decode(changed),
        caption: [Bandstructure of ice Ih computed with MACE-MP-0.],
      ) <fig:phonons-bandstructure-ice-ih-mace-mp-0-zoom>],
  ),
)

==== Timing

#large_box(
  grid(
    columns: 2,
    gutter: 5pt,
    align: horizon,
    figure(
      tablem(
        ignore-second-row: false,
        [
          |supercell|device|optimization time|forces time|
          |1|cuda|45s|8s|
          |2|cuda|44s|50s|
          |3|cuda|44s|21s|
          |4|cpu|5m 23s|1h 47m 7s|
        ],
      ),
      caption: [Execution times with phonopy.],
    ),
    figure(
      tablem(
        ignore-second-row: false,
        [
          |supercell|time|device|
          |2|1m 30s|cuda|
          |4|fail (out of memory)|cuda|
          |4|7h 32m| cpu|
        ],
      ),
      caption: [Execution times with Phonons by @ase.],
    ),
  ),
)

=== Phonons DOS

#large_box(
  grid(
    columns: 2,
    gutter: 5pt,
    figure(
      image("thesis/imgs/holzapfel2021_fig21.png"),
      caption: [
        Reference phonons DOS, taken from @holzapfelCoherentThermodynamicModel2021.
        Comparison of the experimental phonons DOS (green curve) with the theoretical phonons DOS (blue curve) and the neutron scattering function (red curve with roughly adjusted scale).
      ],
    ),
    figure(
      image("simulazioni/02_water/04_crystal_phonons/phonopy/Grafici/dos_mace-ice13-1_s3_mesh=32_zoom.svg"),
      caption: [
        Calculated phonons DOS, using MACE-ICE13-1 and smearing width $sigma=0.05$.
      ],
    ),
  ),
)

As anticipated in the analysis of the bandstructure,
frequencies calculated with MACE-ICE13-1 are shifted toward higher values,
compared to reference data.

#image("simulazioni/02_water/04_crystal_phonons/phonopy/mace_ice13_1_s3_dos.svg")

=== Heat capacity

@flubacherHeatCapacityIce1960 treats the heat capacity of ice at low temperatures.
@holzapfelCoherentThermodynamicModel2021 provides an all-round thermodynamic model of ice Ih, detailing the analysis of heat capacity with 1 Debye and 7 Einstein terms; the analysis of heat capacity considering harmonic terms is reproduced along with reference data in @fig:heat-capacity-mace-holzapfel.

In the quasi-harmonic approximation the lattice vibrations are assumed to be harmonic but with frequencies dependent upon the volume. @leadbetterThermodynamicVibrationalProperties1965

#figure(
  image("simulazioni/02_water/04_crystal_phonons/phonopy/heat_capacity_all_temps.svg"),
  caption: [
    Heat capacity of ice Ih.
    Comparison of results from simulation with different calculators (S3 indicates supercell 3x3x3, otherwise supercell is 2x2x2) and reference data @holzapfelCoherentThermodynamicModel2021.
  ],
) <fig:heat-capacity-mace-holzapfel>

=== Band structure and DOS of $#ce("D2O")$

A further study was performed to analyze the performance of the calculator compared with reference data on deuterated water. @strasslePhononDispersionIce2004
In this scenario, the fidelity is estimated to be worse than the previous case.
The same considerations on the quality of the results hold as above.

#large_box(
  grid(
    columns: (1fr, 0.95fr),
    gutter: 15pt,
    [
      #figure(
        [
          #image("simulazioni/02_water/04_crystal_phonons/phonopy/MACE-ICE13-1/D2O-Ih/band_structure.svg")
          #image("simulazioni/02_water/04_crystal_phonons/phonopy/MACE-ICE13-1/D2O-Ih/total_dos.svg")
        ],
        caption: [
          Band structure and DOS calculated with MACE-ICE13-1.
        ],
      )
    ],
    [
      #figure(
        [
          #image("simulazioni/02_water/04_crystal_phonons/phonopy/MACE-ICE13-1/D2O-Ih/bandstructure_strassle.png")
          #image("simulazioni/02_water/04_crystal_phonons/phonopy/MACE-ICE13-1/D2O-Ih/dos_strassle.png")
        ],
        caption: [
          Band structure and DOS from reference @strasslePhononDispersionIce2004.
        ],
      )],
  ),
)

== MD

=== RDF

Constant NVT @md simulations with Langevin thermostat #footnote[https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.langevin
] were performed under varying external conditions.
A thermostat couples the system to an external heath bath.
The @rdf of the thermalized states is shown in @fig:rdf, compared with reference data @skinnerBenchmarkOxygenoxygenPairdistribution2013 from X-ray diffraction experiment.
The simulated phisical time shall not be less than 100ps, to allow recombination of bonds in the liquid.
Constant NPT simulations should be more appropriate for the computation of physical properties, but they are missing at the present time.
Further analysis can also be made on the study of the diffusion coefficient and the density of the system.

#figure(
  image("simulazioni/02_water/05_md/Grafici/rdf_oo_mace-ice13-1_100ps_nbins=40.svg"),
  caption: [Radial distribution function of oxygens in liquid water.],
) <fig:rdf>

#large_figure(
  grid(
    columns: 2,
    image("simulazioni/02_water/05_md/Grafici/rdf_oo_mace-mp-0_NVT_T=297.15_t=5ps.svg"),
    image("simulazioni/02_water/05_md/Grafici/rdf_oo_mace-ice13-1_nbins=40.svg"),
  ),
  caption: [Comparison of the #glspl("rdf") obtained from @md simulations of liquid water using MACE-MP-0 and MACE-ICE13-1.],
)

#figure(
  image("simulazioni/02_water/05_md/Grafici/temperature_NVT.png"),
)

= Tools
Ibisco, MACE @Batatia2022mace @Batatia2022Design, ASE.

== ASE

Introducing @ase.
The `Atoms` object contains the positions of the atoms and the properties of the cell.

== MACE

=== Fine-tuning a custom model

// TODO Cita kaur2024

Fine-tuning a MACE model is composed of three steps.
// TODO Spiega meglio cosa è e a cosa serve il fine-tuning

1. *Sample the phase space* to obtain representatives of the system under the thermodynamic conditions we are interesed in.
  This was obtained employing NPT dynamics as provided by @ase to generate 10000 images with variety.
  Of this total, 50 representatives were randomly chosen from the images after thermalization.
  See @fig-finetune-sample-temperature for the thermalization pattern.
  // TODO menziona scelta tra termostato di Berensden e Parrinello, con citazioni.
2. *Compute the reference* values for the energies, forces and stresses.
  For this step we executed single point self-consistent calculations on each representative of the samples,
  using a @paw @pbe @dft pseudo-potential through @vasp.
  The output files are then converted, joined and shuffled to compose the training and test sets for the machine learning procedure.
3. *Fine tune* the foundation model on the new dataset; in our case MACE-MP-0 was chosen.
  #footnote([
    Reference procedures are available at https://github.com/ACEsuit/mace?tab=readme-ov-file#finetuning-foundation-models
  ])
  The training and test sets are the input data for the training of the neural network that underlies MACE.
  The training spans 2000 epochs, at the end of which the compiled model is obtained, along with statistics.
  Plot showing neural network loss and error on energy versus epochs are available in @fig-finetune-epochs.

#figure(
  [
    #let my-node(..args) = node(
      corner-radius: 5pt,
      ..args,
    )

    #diagram(
      node-stroke: 1pt,
      my-node((0, 0), [Generate samples \ from NPT MD]),
      edge(
    "-|>",
    [data_for_train.extxyz],
    bend: 45deg,
    // label-sep: 1em,
  ),
      my-node((1, 0), [Compute reference \ DFT pseudopotential]),
      edge(
        "-|>",
        align(
          center,
          [training_set.xyz \ test_set.xyz],
        ),
        bend: 45deg,
      ),
      my-node((2, 0), [Fine tune \ MACE]),
      edge("-|>"),
      my-node((3, 0.5), [Model]),
      edge((2, 0), (3, -0.5), "-|>"),
      my-node((3, -0.5), [Errors]),
    )
  ],
  caption: [Flow diagram for fine-tuning of MACE.],
)

#figure(
  image("tutorial-fine-tuning/1.generate-training/run_2024-06-08/Figure_1_i_T.png"),
  caption: [
    Temperature during sample generation.
  ],
) <fig-finetune-sample-temperature>

#large_figure(
  grid(
    columns: 2,
    image("tutorial-fine-tuning/analysis/loss_over_epochs.svg"),
    image("tutorial-fine-tuning/analysis/mae_e_per_atom_over_epochs.svg"),
  ),
  caption: [
    Loss (left) and @mae of the energy (right) over epochs plots for the fine-tuning of a new model on ice Ih, with MACE-MP-0 small as foundation model.
  ],
) <fig-finetune-epochs>

= Conclusions

The use of #glspl("mlp") allows rapid prototyping through a fast and efficient simulation of representative structures, with a very convenient accuracy/cost ratio.
The combination of foundation models with the ability to fine-tune models to better reproduce the behaviour of particular structures of one's interest, makes this new approach very versatile, and a precious addition in the toolbox of the material scientists.

#box(
  stroke: blue + 2pt,
  inset: 1mm,
  [The accuracy of the new approaches sets the stage for future studies of kinetic effects as wess as full $p-T$ phase diagrams in a reliable and computationally efficient manner. @kapilCompleteDescriptionThermodynamic2022
  ],
)

= Acknowledgments
I would like to express my gratitude to Flaviano Della Pia for his invaluable assistance in learning the new tools and for providing me with insightful advice and practical tutorials. I am also indebted to Andrea Zen for his meticulous guidance and prompt advice. Finally, I would like to thank Dario Alfè for imparting me with the fundamentals of knowledge in this field and for steering my thesis in the right direction.
// This work has been funded by project code PIR01_00011 “IBISCo”, PON 2014-2020, for all three entities (INFN, UNINA and CNR).

#show heading.where(level: 1): it => [
  #[] <empty-page>
  #pagebreak(to: "even", weak: true)
  #[] <new-chapter>
  #pagebreak(to: "odd", weak: true)
  #set text(font: "Neo Euler")
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

