#import "@preview/whalogen:0.2.0": ce
#import "@preview/glossarium:0.4.1": make-glossary, print-glossary, gls, glspl
#import "@preview/hydra:0.4.0": hydra
#import "@preview/fletcher:0.4.5" as fletcher: diagram, node, edge
#import "@preview/tablem:0.1.0": tablem
#import "@preview/physica:0.9.3": *
#show: make-glossary

// For glossarium links
#show link: set text(fill: blue.darken(60%))

#let title = "Properties of molecular crystals using machine learning potentials"
#set document(title: title, author: "Mariano Mollo")
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

  #text(1.5em, [*UNIVERSITÀ DEGLI STUDI DI NAPOLI \ "FEDERICO II"*])

  #v(3mm)

  // University Logo
  #image("thesis/imgs/University_Federico_II_Logo.svg", width: 25%)

  #v(1cm)

  *Scuola Politecnica e delle Scienze di Base*

  *Area Didattica di Scienze Matematiche Fisiche e Naturali*

  #v(8mm)

  *Dipartimento di Fisica "Ettore Pancini"*

  #v(20mm)

  _Laurea Magistrale Sperimentale in Fisica_

  #v(5mm)

  #text(1.5em, title)

  #v(25mm)

  #grid(
    columns: 2,
    align: (left, right),
    column-gutter: 1fr,
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

#par(justify: true, first-line-indent: 0pt)[
  = Abstract

  Molecular crystals play an important role in the field of materials science,
  particularly in drug development, electronics, and renewable energy sectors.

  In this work, we will use recently developed #glspl("mlp") to model the structure and dynamics of molecular crystals along with their thermodynamic stability, using water as a showcase system.
  Water is ubiquitous in nature and of fundamental relevance to physics, biology, geology, materials science and engineering.
  Its numerous anomalies, arising from the delicate interplay of hydrogen bonding and dispersion forces, make it a hard test for computational approaches.

  Traditional approaches often grapple with the trade-off between computational
  expense and accuracy. The application of #glspl("mlp") captures
  complex intermolecular interactions with the accuracy of ab initio approaches but at a much cheaper computational cost.

  In this work, we will test different @mlp models for the prediction of properties of molecular crystals, including lattice energy, phonons dispersion, and finite temperature dynamics.
]

#pagebreak()

#[#show outline.entry.where(level: 1): it => {
    v(12pt, weak: true)
    strong(it)
  }
  #outline()
]

#[
  // Vorrei troncare le entry troppo lunghe
  #show outline.entry: it => {
    it
  }
  #outline(
    title: [List of Figures],
    target: figure.where(kind: image),
  )
]

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
    (key: "ml", short: "ML", long: "Machine Learning"),
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
    (
      key: "hdnnp",
      short: "HDNNP",
      long: "high-dimensional neural network potential",
    ),
    (key: "mpnn", short: "MPNN", long: "Message Passing Neural Network"),
    (key: "gnn", short: "GNN", long: "Graph Neural Network"),
    (
      key: "ibisco",
      short: "IBiSCo",
      long: "Infrastructure for BIg data and Scientific Computing",
    ),
    (key: "cnn", short: "CNN", long: "Convolutional Neural Network"),
    (key: "adaline", short: "Adaline", long: "ADAptive LInear Neuron"),
    (key: "nn", short: "NN", long: "Neural Network"),
    (key: "ks", short: "KS", long: "Kohn-Sham"),
    (key: "lda", short: "LDA", long: "Local Density Approximation"),
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
#show heading.where(level: 3): it => [
  #it
  #v(4mm)
]
#show heading.where(level: 4): it => [
  #it
  #v(3mm)
]

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

Molecular crystals represent a significant area of study within materials science due to their diverse applications in fields such as pharmaceuticals, electronics, and renewable energy.
Computational approaches play an important role in molecular crystals research, as their application can aid experiments toward prediction of stable phases with targeted properties.
In drug development, for instance, the properties of molecular crystals directly influence the effectiveness and delivery of medications, as exemplified in @fig:paracetamol-tableting. @blagdenCrystalEngineeringActive2007
Beyond pharmaceuticals, these materials are key to advancing technologies involving semiconductors and energy storage systems, given their tunable electronic and optical properties. @reddyMechanicalPropertiesMolecular2010 @martinsTemperaturePressureInducedProton2009 @karkiImprovingMechanicalProperties2009 @reillyUnderstandingRoleVibrations2013

#figure(
  image("thesis/imgs/karki2009_graphical_abstact.jpg", width: 50%),
  caption: [
    Poor mechanical properties of paracetamol are improved through the strategy of cocrystal formation.
    Graphical Abstract taken from @karkiImprovingMechanicalProperties2009.
  ],
) <fig:paracetamol-tableting>

Despite the critical importance @priceComputationalPredictionPharmaceutical2004 of understanding and predicting the properties of molecular crystals, traditional computational approaches face challenges in balancing accuracy with computational cost.
Methods such as @dft offer accurate insights into molecular interactions but are prohibitively expensive for large-scale simulations.
Conversely, classical force field methods scale well but often lack the precision needed for complex systems.
This trade-off has driven the exploration of alternative methods that combine the best of both worlds.

#figure(
  image("thesis/imgs/images_large_cr0c00868_0001.jpeg", width: 80%),
  caption: [
    From @behlerFourGenerationsHighDimensional2021[Figure 1].
    Pyramid of potentials for atomistic simulations illustrated for the example of water and aqueous systems.
    Using high-level wave function-based methods as represented by Ψ only the geometries of small systems such as water clusters in vacuum are accessible, while DFT is the standard method to determine simple properties such as RDFs of liquid water in ab initio molecular dynamics simulations.
    Very large-scale simulations of complex systems such as electrolytes or solid−liquid interfaces, or the determination of complex thermodynamic properties, can only be carried out using atomistic potentials such as forceﬁelds.
    Also MLPs allow the study of these systems.
  ],
)

Recently, #glspl("mlp") have emerged as a promising solution.
These approaches leverage on advanced neural network architectures to model the potential energy surfaces of molecular systems with a level of accuracy comparable to ab initio methods but at a fraction of the computational cost.
By training on existing quantum mechanical data, #glspl("mlp") can predict intermolecular interactions and dynamic behaviours with high fidelity, making them suitable for studying complex systems like molecular crystals.
@batatiaFoundationModelAtomistic2023
@schranMachineLearningPotentials2021
@bjorneholmWaterInterfaces2016
@kapilCompleteDescriptionThermodynamic2022

Since their introduction about 25 years ago, #glspl("mlp") have become an important tool in the ﬁeld of atomistic simulations.
After the initial decade, in which neural networks were successfully used to construct potentials for rather small molecular systems, the development of #glspl("hdnnp") in 2007 opened the way for the application of ML potentials in simulations of large systems containing thousands of atoms.
@behlerFourGenerationsHighDimensional2021
Depending on the specific task it has been estimated that mixed ab initio and @mlp calculations require between three and ten times less the core hours of purely ab initio methods @kapilCompleteDescriptionThermodynamic2022[p. 5], while state of the art #glspl("mpnn", long: true) achieve speedups of up to $10^5$ times less the compute time on the prediction of physical properties of molecules (see @fig:mpnn-speedup).

#figure(
  image("thesis/imgs/gilmerNeuralMessagePassing2017_Figure1.png"),
  caption: [
    A Message Passing Neural Network predicts quantum properties of an organic molecule by modeling a computationally expensive DFT calculation.
    Image taken from @gilmerNeuralMessagePassing2017.
  ],
) <fig:mpnn-speedup>

This thesis investigates the capabilities of #glspl("mlp"), with a specific focus on modeling the properties of water as a molecular crystal.
Water is ubiquitous in nature and, as stated above, of fundamental relevance to physics, biology, geology, materials science, and engineering, that serves as an ideal test case due to its well-documented polymorphic behaviours and the wealth of literature available for benchmarking. @schranMachineLearningPotentials2021 @chengInitioThermodynamicsLiquid2019
Studying its numerous anomalies---arising from the delicate interplay of hydrogen bonding and dispersion forces---is a challenge for computational approaches, and allows us to discover new physics and advance various scientific applications.
The research primarily explores the accuracy of #glspl("mlp") in predicting molecular structural stabilities, lattice energies of different ice polymorphs, and dynamical properties of crystals, such as phonon spectra.

// TODO kapil 29 : V. Kapil, E. Engel, M. Rossi, M. Ceriotti, Assessment of approximate methods for anharmonic free energies. J. Chem. Theory Comput. 15, 5845–5857 (2019).

Despite significant advancements, state-of-the-art #glspl("mlp") still suffer from several limitations, such as the large amount of data needed for training and the limited number of chemical species that can be included in the model.
However, extremely recent @mlp architectures, such as MACE, have made overcoming these limitations possible, leading to the development of foundational models for chemistry and materials science. @batatiaFoundationModelAtomistic2023 @dengCHGNetPretrainedUniversal2023 @liKohnShamEquationsRegularizer2021 @cheonDatasetRandomRelaxations2023

The following chapters will introduce the theoretical foundations and the tools needed to pursue this question.
@sec:theory briefly describes the physical priciples and the algorithms on which #glspl("mlp") like MACE are built onto.
@sec:results-1 and @sec:results-2 show the results obtained through the test of the new MACE calculator on known atomic configurations and discuss its performances.
@sec:tools details the tools employed to run the computer simulation experiments;
the algorithms section, @sec:mace, mirrors and details the practical implementation of the objects first outlined in @sec:gnn.

= Theory <sec:theory>

In this chapter, we explore the theoretical background and methodologies that underpin the analysis of molecular crystals using #glspl("mlp", long: true).
The study of molecular crystals involves various computational and physical models that allow us to predict their structure, dynamics, and thermodynamic stability.
The chapter is orgnanized into several key areas:

#let threed-harmonic-solids-title = "Three dimensional harmonic crystalline solids"
#let dft-title = "Density Functional Theory"
#let md-title = "Molecular Dynamics"
#let thermodynamics-crystal-stability-title = "Thermodynamics of crystal stability"
#let machine-learning-title = "Machine Learning"
#let gnn-title = "Graph Neural Networks"

- #strong(threed-harmonic-solids-title), @sec:3d-harmonic-solids:
  We introduce the harmonic approximation for crystalline solids and explain its relevance to modelling vibrational properties through normal modes, focusing on phonons and their contribution to the Helmholtz energy.
- #strong(dft-title), @sec:dft:
  @dft serves as the core quantum mechanical approach used to compute the electronic properties and structural optimizations of materials.
  Here, we describe its formulation and the Kohn-Sham method for solving the many-body Schrödinger equation.
- #strong(md-title), @sec:md:
  This section covers the simulation of atomic and molecular motion using classical mechanics.
  We explain the Verlet algorithm, thermostats for temperature control, and how these methods help model the dynamic behaviour of molecular systems at finite temperatures.
- #strong(thermodynamics-crystal-stability-title), @sec:thermodynamics-crystal-stability:
  We detail the calculation of lattice energy, which measures the stability of crystals, and discuss the important role of dispersion interactions, including methods like DFT-D3, in capturing non-covalent forces in crystalline systems.
- #strong(machine-learning-title), @sec:machine-learning:
  Machine Learning techniques, particularly neural networks, are introduced as modern tools for predicting molecular properties with reduced computational cost.
  We explore single-layer and multi-layer neural network architectures and their training processes.
- #strong(gnn-title), @sec:gnn:
  As an advanced form of neural networks, #glspl("gnn") are discussed, with emphasis on how they enable message-passing architectures for modelling molecular systems.
  We describe how #glspl("gnn") can be employed to capture atomic interactions in a computationally efficient manner.

Each section provides the necessary theoretical foundation for understanding the subsequent results and discussions, with a focus on applying these methods to study the properties of molecular crystals and, specifically, ice polymorphs.
The content found from @sec:3d-harmonic-solids to @sec:md, is significantly based upon the lecture notes on Statistical and Computational Physics, by D. Alfè, studied during my masters degree. @alfeNotesStatisticalComputational2023

== #threed-harmonic-solids-title <sec:3d-harmonic-solids>
The crystal can be described by a collection of independent harmonic oscillators with potential energy:
$
  U = U_0 + sum_(i=1)^(3N) 1 / 2 M omega_i^2 q_i^2,
$ <eq:termcomp-7.1>
where $M$ is the mass of the particles, $q_i$ a set of _normal coordinates_, and $U_0$ the energy of the system in its ground state.
This approximation is known as the _harmonic approximation_.
The Newton's equations of motion for the normal coordinates are:
$
  M dv(q_i, t, 2) = - pdv(U, q_i) = - M omega_i^2 q_i.
$

Consider a three dimensional system made of particles arranged on a _periodic lattice_, at equilibrium positions ${va(r^0)} := va(r_1^0), dots, va(r_N^0), dots$.
This _Bravais lattice_ is completely defined by a set of three _lattice vectors_, $va(a_1), va(a_2), va(a_3)$, and every point on the lattice can be obtained as:
$
  va(r_j^0) = n va(a_1) + m va(a_2) + l va(a_3),
$
where $j$ is a shorthand for the three integers $n, m, l$.
We also define the _reciprocal lattice vectors_:
$
  va(b_1) := 2 pi (va(a_2) times va(a_3)) / V, quad
  va(b_2) := 2 pi (va(a_3) times va(a_1)) / V, quad
  va(b_3) := 2 pi (va(a_1) times va(a_2)) / V,
$
where $V := (va(a_1) times va(a_2)) dot va(a_3)$ is the volume of the _primitive cell_.
This is the smallest unit of volume that can be periodically repeated to fill the whole space, and clearly it can be constructed in multiple (infinite) ways.
A convenient construction is the _Wigner-Seitz_ primitive cell, which is the locus of points in space closer to a particular lattice point than to any other lattice point;
this procedure is the application of Voronoi decomposition to a crystal lattice.
The reciprocal lattice vectors span the _reciprocal space_, which is also a periodic lattice with points at positions:
$
  va(g) := n' va(b_1) + m' va(b_2) + l' va(b_3),
$
with $n', m', l'$ integer numbers.
The Wigner-Seitz cell in reciprocal space is called the _first Brillouin zone_.
The reciprocal vectors and the lattice vectors satisfy the orthogonality relations:
$
  va(b_i) dot va(a_j) = 2 pi delta_(i j).
$
Each lattice site does not need to be occupied by just one particle; indeed, the generalization to crystals with more than one particle per lattice site is straightforward.

Let $U_({va(r^0)}) ({va(r)})$ be the potential energy function, defined by the interacting particles when they are in their equilibrium positions ${va(r^0)}$.
With this we mean that if we displace the particles by amounts $va(u_i) := va(r_i) - va(r_i^0)$, we assume that these displacements are not changing the functional form of the potential, but only the computed value changes.
This assumption can only be valid if the displacements ${va(u)}$ are small.
An extreme case where this cannot possibly hold is when the atoms move so much that the crystal assumes a different crystal structure, or it melts.
For zero displacements, when all particles are in their equilibrium positions, the potential energy can be computed from knowledge of the lattice vectors only, i.e., one only needs information about the primitive cell.
If instead the particles are displaced from their equilibrium positions, then we need the positions of all of them.
For small enough displacements ${va(u)}$ the potential energy can be expanded around its minimum, where all particles are in their equilibrium positions ${va(r^0)}$:
$
  U({va(r)})
  = U_0
  + 1 / 2 sum_(i,j) va(u_i) dot Phi (
    va(r_i^0) - va(r_j^0)
  ) dot va(u_j) + cal(O)(u^3),
$
where $U_0 := U({va(r^0)})$, and the linear term is absent because we are expanding around the minimum of the potential.
The _force constants matrix_ $Phi$ is defined as:
$
  Phi (va(r_i^0) - va(r_j^0)) := (
    pdv(U({va(r)}), va(r_i), va(r_j))
  )_(va(r_i) = va(r_i^0) \ va(r_j) = va(r_j^0))
  = mat(
    phi_(i j)^(x x), phi_(i j)^(x y), phi_(i j)^(x z);
    phi_(i j)^(y x), phi_(i j)^(y y), phi_(i j)^(y z);
    phi_(i j)^(z x), phi_(i j)^(z y), phi_(i j)^(z z);
  ) \
  phi_(i j)^(alpha beta) = pdv(U({va(r)}), alpha_i, beta_i),
$ <eq:force-constants-matrix>
where $alpha_i$ and $beta_j$ run over the three cartesian components of $va(r_i^0)$ and $va(r_j^0)$.
If the displacements are not too large and the $cal(O)(u^3)$ terms can be ignored, the force acting on particle at position $va(r_i^0)$ due to the displacements $va(u_j)$ of all particles in the system, including $va(u_i)$, is:
$
  va(f_i) = -sum_j Phi(va(r_i^0) - va(r_j^0)) dot va(u_j).
$

The dependence of the force constants matrix on the difference $va(r_i^0) - va(r_j^0)$ rather than on $va(r_i^0)$ and $va(r_j^0)$ separately is due to translational invariance.
The force constants matrix satisfies other important properties:
+ $Phi(va(r_i^0) - va(r_j^0)) = Phi(va(r_j^0) - va(r_i^0))$, which is implied by the fact that the double derivative in @eq:force-constants-matrix is invariant upon changing the order of differentiation.
+ $sum_j Phi(va(r_j^0)) = 0$. This property can be understood by imagining to displace the whole crystal rigidly by some vector $va(c)$.
  Clearly, such a displacement cannot affect the forces acting on the particles, because the relative differences between them are unaffected by the translation, and so we must have $va(f_i) = -sum_j Phi(va(r_i^0) - va(r_j^0)) dot va(u_j) = -sum_j Phi(va(r_i^0) - va(r_j^0)) dot (va(u_j) + va(c))$.
  Therefore, $sum_j Phi(va(r_i^0) - va(r_j^0)) dot va(c) = 0$.
  Since this is independent on the particular choice of $va(r_i^0)$ and $va(c)$, we must have $sum_j Phi(va(r_j^0)) = 0$.

We can obtain the normal modes of the system and show that the potential energy is given by the sum of squares of the normal modes.
To obtain the dispersion relation, consider the Newton's equation of motion for a particle at one particular lattice position, $va(r_i^0)$:
$
  M dv(va(u_i), t, 2) = va(f_i).
$ <eq:termcomp-7.54>
The motion of the particle around its equilibrium position can be described by the form:
$
  va(u_i)(t) prop va(epsilon) e^(i (va(q) dot va(r_i^0) - omega t)),
$ <eq:termcomp-7.55>
where $va(epsilon)$ is a _polarization vector_, which defines the direction of oscillation of the particle, and $va(q)$ is a _wave vector_.
The physical motion of the particles is obtained by taking the real part of @eq:termcomp-7.55.
Substituting @eq:termcomp-7.55 into @eq:termcomp-7.54 we obtain:
$
  M omega^2 va(epsilon)
  = sum_j e^(i va(q) dot (
    va(r_j^0) - va(r_i^0)
  )) Phi(va(r_i^0) - va(r_j^0)) dot va(epsilon).
$ <eq:termcomp-7.56>
The sum over $j$ runs over all lattice sites, and so it is independent on our choice of $va(r_i^0)$.
We can therefore choose $va(r_i^0) = va(0)$ and replace the difference $va(r_j^0) - va(r_i^0)$ simply with $va(r_j^0)$.

=== Phonons
If we introduce
the _dynamical matrix_:
$
  D(arrow(q)) :=
  1 / M sum_j e^(i arrow(q) dot arrow(r)_j^0) Phi(arrow(r)_j^0)
$ <eq:dynamical-matrix>
we can rewrite @eq:termcomp-7.56 as:
$
  omega^2 arrow(epsilon) = D(arrow(q)) dot arrow(epsilon)
$ <eq:dynamical-matrix-eigenvalue-problem>
The transformation in @eq:dynamical-matrix is also known as a _lattice Fourier transform_.
Note that it is sufficient to define the dynamical matrix in the first Brillouin zone, because any vector in reciprocal space can be written as the sum of a vector in the Brillouin zone plus a reciprocal lattice vector.
We have:
$
  e^(i (va(q) + va(g)) dot va(r_j^0))
  = e^(i va(q) dot va(r_j^0)) e^(i va(g) dot va(r_j^0))
  = e^(i va(q) dot va(r_j^0)),
$
because:
$
  va(g) dot va(r_j^0) = (n n' + m m' + l l') 2 pi,
$
and so the dynamical matrix has the same periodicity of the reciprocal lattice.

The property $sum_j Phi(va(r_j^0)) = 0$ implies $D(va(q)) = va(0)$; therefore, the dispersion relation gives zero frequencies at zero wavevector.
This _sum rule_ is very useful as a practical sanity check of the consistency of the calculations.

To solve @eq:dynamical-matrix-eigenvalue-problem we need to solve an eigenvalue problem.
We need to find a reference frame in which the matrix $D(va(q))$ is diagonal.
The elements of the diagonal, $omega_(va(q), s)^2$, with $s in {1,2,3}$ representing the _branch number_, are the _eigenvalues_, and the three vectors $va(epsilon_(va(q), s))$, the polarizations of the normal modes, are the _eigenvectors_ that define the reference frame.

Each normal mode represents a collective oscillation of the particles in the system with frequency $omega_(va(q), s)$.
These collective oscillations are called _phonons_.

The potential energy can be written in the form in @eq:termcomp-7.1 and so the partition function has the form of either @eq:termcomp-7.7 or @eq:termcomp-7.10, depending on if the system is treated classically or quantum-mechanically:
$
  Z = e^(-beta U_0) product_(i=1)^(3N) (k_B T) / (planck.reduce omega_i), quad "in the classical limit",
$ <eq:termcomp-7.7>
and
$
  Z = e^(-beta U_0) product_(i=1)^(3N) Z_i, quad "in the quantum description",
$ <eq:termcomp-7.10>

with each term in the product having the form:
$
  Z_i = sum_(n=0)^infinity e^(-E_n^i / (k_B T))
  = sum_(n=0)^infinity e^(- (n + 1 / 2) (planck.reduce omega_i) / (k_B T))
  = (e^(-1 / 2 (planck.reduce omega_i) / (k_B T))) / (1 - e^(- (planck.reduce omega_i) / (k_B T))).
$

=== The Helmholtz energy
// Da §4 relazione termodinamica computazionale

The harmonic contribution of phonons to the Helmholtz energy per atom is represented by the equation:

$
  &F_"harm" (V,T) = \
  =&
  1 / Omega integral_Omega sum_(s=1)^3 (planck.reduce omega_(arrow(q),s)) / 2 dif q
  + k_B T 1 / Omega integral_Omega sum_(s=1)^3 ln(1-exp(-(planck.reduce omega_(arrow(q),s))/(k_B T))) dif q
$

In the high temperature limit, defined by $(planck.reduce omega_(arrow(omega), s))/(k_B T) << 1$, the harmonic term reduces to the classical expression:

$
  F_"harm"^c = (k_B T) / Omega integral_Omega sum_(s=1)^3 ln((planck.reduce omega_(arrow(q),s))/(k_B T)) dif q
$

For computation efficiency purposes, the integrals above are approximated with a sum running on points sampled in the Brillouin zone, denoted by BZ in the sum:

$
  1 / Omega integral dif q approx 1 / (N_(arrow(q))) sum_(arrow(q) in "BZ")
$

Using this approximation, the Helmholtz energy per atom becomes:

$
  F (V,T)
  &= U_0(V) \
  &+ 1 / (N_va(q)) sum_(s=1)^3 sum_(va(q)) [
    (hbar omega_(va(q),s) (V)) / 2
    + k_B T ln(
      1 - exp(
        - (hbar omega_(va(q), s) (V)) / (k_B T)
      )
    )
  ],
$ <eq:termcomp-7.63>
with $U_0$ the energy per atom of the system in its ground state.
In the classical limit, the expression becomes:

$
  F_"classical" (V,T) = U_0(
    V
  ) + (k_B T) / (N_va(q)) sum_(s=1)^3 sum_(arrow(q)) ln (
  planck.reduce omega_(arrow(q),s)(V)) / ( k_B T
  )
$ <eq:termcomp-7.64>

The sum over $va(q)$ runs over all vectors in the Brillouin zone, which are infinite for a Bravais lattice that contains an infinite number of sites.
Indeed, the Helmholtz energy of a crystal with an infinite number of particles would also be infinite.
The physical quantity of interest is indeed the Helmholtz energy per particle; this is the reason why the sums in @eq:termcomp-7.63 and @eq:termcomp-7.64, using a finite grid of $N_(va(q))$ points in the Brillouin zone and dividing the total Helmholtz energy by $N_(va(q))$, yield finite values.
Both the ground state energy and the frequencies explicitly depend on the volume $V$.
The relationship linking the frequencies with volume and temperature is responsible for the different temperature dependence of the Helmholtz energy at different volumes; these terms cause the phenomenon of thermal expansion in solids.

=== The small displacement method <sec:small-displacement-method>
// Vedi anche §3.4.2 Della Pia, §7.1.3 Alfè

Let us consider the crystal in the ground state and displace one particle from its equilibrium position.
Let the displacement be along the $x$ axis, $arrow(u) := (u, 0, 0)$.
The forces acting on all particles in the system, including the displaced one, are given by:

$
  (f_i^x, f_i^y, f_i^z)
  = -Phi(arrow(r)_i^0) dot (u, 0, 0)
  = -u vec(phi_(0i)^(x x), phi_(0i)^(y x), phi_(0i)^(z x))
$

from which we obtain the first columnt of the force constant matrix:

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
when this happens, the only term contributing appreciably in the sum in @eq:force-constants-matrix-supercell is that with $arrow(L) = arrow(0)$;
as we increase the size of the supercell, the supercell force constants matrix asymptotically approaches the force constants matrix,
$Phi_"SC" (arrow(r)_i^0) tilde.eq Phi(arrow(r)_i^0)$.
In this limit, the Fourier interpolation is accurate everywhere.

== #dft-title <sec:dft>
#gls("dft", long: true) is in principle an exact formulation of the many-body electron problem, in terms of the ground state density of the electrons rather than the ground state wavefunction.
One of its numerous applications is the determination of structural and electronic properties of materials in solid-state physics.
In this thesis, the @dft calculations are performed with the #gls("vasp", long: true). @kresseInitioMolecularDynamics1993 @kresseEfficiencyInitioTotal1996 @kresseEfficientIterativeSchemes1996
Here, we briefly describe the underlying theory of @dft.

Given the wavefunction $Psi (va(r_1), dots, va(r_N))$ and the Hamiltonian $hat(H)$ of the time-independent many-body Schrödinger equation for a system of $N$ electrons, not including spin variables, $hat(H) Psi (va(r_1), dots, va(r_N)) = E Psi (va(r_1), dots, va(r_N))$, the Hamiltonian is given by the sum of the kinetic, external potential, and electron-electron operators, $hat(H) = hat(T) + hat(V)_"ext" + hat(V)_"ee"$, defined by:
$
  hat(T) Psi (va(r_1), dots, va(r_N))
  := - 1 / 2 sum_(i=1)^N laplacian_i Psi (va(r_1), dots, va(r_N)),
$
$
  hat(V)_"ext" Psi (va(r_1), dots, va(r_N))
  := sum_(i=1)^N v(va(r_i)) Psi (va(r_1), dots, va(r_N)),
$
with $v(va(r_i))$ the value of the external potential felt by electron $i$ at position $va(r_i)$, and
$
  hat(V)_"ee" Psi (va(r_1), dots, va(r_N))
  := sum_(i,j=1 \ i<j) 1 / (|va(r_i) - va(r_j)|) Psi (va(r_1), dots, va(r_N)).
$

If the wavefunction $Psi$ is normalized over the volume of the system $V$,
$
  integral_V dif^3 va(r_1) dots dif^3 va(r_N)
  Psi^star (va(r_1), dots, va(r_N)) Psi (va(r_1), dots, va(r_N))
  = 1,
$
then the energy $E$ is obtained from:
$
  E =
  integral_V dif^3 va(r_1) dots dif^3 va(r_N) Psi^star (
    va(r_1), dots, va(r_N)
  ) hat(H) Psi (va(r_1), dots, va(r_N))
$

Computing the electron density gives:
$
  rho(va(r)) = N integral_V dif^3 va(r_2) dots dif^3 va(r_N) |Psi(va(r), va(r_2), dots, va(r_N))|^2.
$ <eq:electron-density>

The total potential energy due to the external potential is a _functional_ of the electron density, which we indicate with the notation $V_"ext" [rho]$:
$
  V_"ext" &=
  integral_V dif^3 va(r_1) dots dif^3 va(r_N) Psi^star (va(r_1), dots, va(r_N))
  sum_(i=1)^N v(va(r_i))
  Psi (va(r_1), dots, va(r_N)) \
  &= sum_(i=1)^N integral_V dif^3 va(r_1) dots dif^3 va(r_N) Psi^star (
    va(r_1), dots, va(r_N)
  ) v(va(r_i)) Psi (va(r_1), dots, va(r_N)).
$

Since the electrons are all identical, these integrals are also all identical, each one equal to:
$
  integral_V dif^3 va(r) dif^3 va(r_2) dots dif^3 va(r_N)
  Psi^star (va(r_1), dots, va(r_N)) v(va(r)) Psi (va(r_1), dots, va(r_N)).
$

Since we have $N$ of them, we obtain:
// Alfè, eq. (8.11)
$
  V_"ext"
  = N integral_V dif^3 va(r) dif^3 va(r_2) dots dif^3 va(r_N)
  Psi^star (va(r_1), dots, va(r_N)) v(va(r)) Psi (va(r_1), dots, va(r_N))
$

This shows, comparing with @eq:electron-density, that the potential energy due to the external potential is a _functional_ of the electron density:
$
  V_"ext" [rho] = integral_V dif^3 va(r) rho(va(r)) v(va(r)).
$

=== The Hohenberg-Kohn theorems

In the Schrödinger approach, both $Psi$ and $rho_0 (va(r))$ are functionals of $v(va(r))$.
The *Hohenberg-Kohn theorems* @hohenbergInhomogeneousElectronGas1964 state that $V_"ext" [rho_0]$, where $rho_0$ is the ground state density, and $Psi_0$, the ground state wavefunction of the system, and thus all the physical properties of the system, are _unique_ functionals of $rho_0$.

In particular, without giving formal proof here, we say that the ground state energy, $E_0$, is a functional of $rho_0$, as well as the kinetic and electron-electron interaction contributions:
$ // Alfè, eq. (8.17)
  E_0 = E[rho_0] = T[rho_0] + V_"ee" [rho_0] + V_"ext" [rho_0].
$
The functional $ F[rho] := T[rho] + V_"ee" [rho] $ does not depend on the external potential, and therefore it is a _universal_ functional of the density.

=== The Kohn-Sham method

Unfortunately, the explicit dependence of the $F[rho]$ functional with respect to $rho(va(r))$ is unknown;
so, an exact solution is not possible.

Kohn and Sham extracted from the universal functional $F[rho]$ the classical Coulombian energy, defining the functional:
$ // Della Pia, eq. 3.32
  G[rho] := F[rho] - 1/2 integral (rho(va(r)) rho(va(r')))/(|va(r) - va(r')|) dif^3 va(r) dif^3 va(r'),
$
where $rho(va(r))$ is a generic electronic density.
The functional $G[rho]$ contains, according to the previous considerations, the quantum contributions of the Coulomb interaction and the kinetic energy of the interacting system of electrons.
Subsequently, for comparison with non-interacting systems of electrons, they defined the exchange-correlation energy functional as:
$ // Della Pia, eq. 3.33
  E_"xc" [rho] := G[rho] - T_s [rho],
$
where $T_s [rho]$ is the kinetic energy functional of the unique non-interacting system of electrons having the same ground state electronic density as the system under consideration.
Hence, $E_"xc" [rho]$ contains the quantum contributions of the Coulomb interaction and the remaining contribution of the kinetic energy of the interacting electrons.

The *Kohn-Sham method* @kohnSelfConsistentEquationsIncluding1965 provides a working procedure to find the ground state density, defining the Kohn-Sham potential, $v_"KS" [rho] (va(r)) := v(va(r)) + integral_V (rho(va(r'))) / (|va(r) - va(r')|) dif^3 va(r') + (delta E_"xc" [rho]) / (delta rho)$, that leads to the Kohn-Sham self-consistent equations:
$ // Eq. (8.39) Alfè, (3.39) Della Pia
  hat(h)_"KS" [rho] (va(r)) psi_n (va(r)) =
  [-1 / 2 nabla^2 + v_"KS" [rho] (arrow(r))] psi_n (arrow(r))
  = epsilon_n psi_n (arrow(r))
$ <eq:kohn-sham-equations>

The equations are coupled through the effective potential $v_"KS" [rho]$, as $rho(va(r))$ depends on all the $psi_n (va(r))$.

Since the orbitals $psi_n$ are ortho-normal, the ground state electron density of the system is obtained from the solution of @eq:kohn-sham-equations[Equations] as:
$
  // Eq. (3.41) Della Pia, (8.23) Alfè
  // logseq://graph/softseq?block-id=66f33133-67b5-46ea-abe5-ee9003827258
  rho_0 (va(r)) = sum_(n=1)^N |psi_n (va(r))|^2.
$ <eq:ground-state-density>

// logseq://graph/softseq?block-id=66f33215-091e-46ec-bd3c-1ce9a47d0558
It is useful to recast the variation of the density in terms of variations of the single particle wavefunctions $psi_n$.
Such a variation has to be performed while keeping the wavefunctions ortho-normal, which gives a set of $M^2$ constraints:
$ // Termcomp eq. 8.28
  integral_V dif^3 va(r) psi_i^star (va(r)) psi_j (va(r)) = delta_(i j).
$
To do that, we define the functional:
$ // Termcomp eq. 8.29
  Omega [ { psi_n }, {epsilon_(i j)} ]
  := E[rho]
  - sum_(i,j=1)^(N/2) 2 epsilon_(i j)
  ( integral_V dif^3 va(r) psi_i^star (va(r)) psi_j(va(r)) - delta_(i j) ),
$
where the terms $epsilon_(i j)$ are the Lagrange multipliers associated to the constraints, and we impose the condition:
$ // Termcomp eq. 8.30
  delta Omega [{psi_n}, {epsilon_(i j)}] = 0.
$

// logseq://graph/softseq?block-id=66f328eb-1dc2-4e46-85f3-68f32ad31ce7
The @ks equations show that the extremes of the functional $Omega$ are obtained for any ensemble of $N/2$ eigenstates of $hat(h)_"KS" [rho]$ ($N$ eigenstates, if we ignore the spin degeneracy), and the ground state is obtained by ﬁnding its minimum of the total energy with respect to any choice of $N/2$ state.
In practice, to approximate the interacting system one always takes the lowest $N/2$ @ks states, although this may not necessarily be the correct choice (see @sec:v-representability).

One can introduce a set of $L$ basis functions, ${phi_m}_(m = 1, dots, L)$, to construct the matrix $epsilon_(i j) = braket(phi_i, hat(h)_"KS" [rho], phi_j)$, and diagonalize it.
The eigenvectors are of the type:
$ // Eq. (8.40) Alfè, (3.42) Della Pia
  psi_n = sum_(i=1)^L c_i^n phi_i
$
If the basis set ${phi_m}$ is complete, the solution is exact.
However, usually this means including an infinite number of elements, which cannot be done in practice;
therefore, the solution to the @ks equations is only approximate.
By increasing the number of basis functions, one can drive the calculations to convergence, where the energy and other properties are obtained within some predefined threshold.
The usual variational principle applies, and so, by including more and more elements in the basis set, the ground state energy decreases monotonically.
The rate of decrease of the energy can be used to judge the level of convergence.

// logseq://graph/softseq?block-id=66f33452-dd0c-4443-86c6-5e64dba07f2c
Since the @ks potential depends on $rho$, @eq:kohn-sham-equations[Equations] have to be solved self-consistently.
This is typically done by iteration, in which one starts with some initial guess for the electron density $rho_1$, constructs $v_"KS" [rho_1]$ and solves @eq:kohn-sham-equations[Equations].
With those solutions construct $rho_2$ using @eq:ground-state-density or more advanced algorithms, and solve @eq:kohn-sham-equations[Equations] again, using the newly constructed $v_"KS" [rho_2]$.
The algorithm runs until the difference between $rho_(j+1)$ and $rho_j$ is below some acceptable threshold.

Multiplying @eq:kohn-sham-equations[Equations] by $psi_n^star (va(r))$, integrating over $va(r)$, summing over $n$, the total energy can be obtained from:
$ // Eq. (8.41) Alfè
  E = 2 sum_(n=1)^(N / 2) epsilon_n - 1 / 2 integral_V (rho (
    arrow(r)'
  ) rho(arrow(r))) / (|arrow(r) - arrow(r)'|) dif^3 arrow(r)' dif^3 arrow(r)
  - integral_V (delta E_"xc") / (delta rho(arrow(r))) rho(arrow(r)) dif^3 arrow(r) + E_"xc" [rho].
$ <eq:termcomp-8.41>

=== v-representability <sec:v-representability>
// Da Lecture Notes di Alfè, p. 150
The condition that guarantees that the first $N/2$ states give the minimum energy in @eq:termcomp-8.41 is known as _v-representability_, that is the ground state density of the interacting system is the same as the ground state density of _some_ non-interacting system.
This also implies that the total energy, @eq:termcomp-8.41, would be exact if one had the exact form of the exchange-correlation functional, $E_"xc"$.
It is possible, however, that the ground state density of a system is not _v-representable_, i.e. there is no non-interacting system whose ground state density is the same as that of the interacting system.
In this case we need to go back to the variational principle, and the $N/2$ states in the non-interacting system that minimize the energy may not be the first $N/2$. @parrDensityFunctionalTheoryAtoms1994

=== The local density approximation

The explicit definition of the exchange-correlation energy is unknown.
We need to approximate the exchange-correlation functional in the Kohn-Sham equations.
The simplest definition of $E_"xc" [rho]$ is based on the assumption that the density in a given point in space, $va(r)$, is a smooth function around this point.
This is strictly true for a homogeneous electron gas, for which the potential $v_"xc" (va(r))$ can be seen to depend on the local density;
however, for an inhomogeneous electron gas it becomes a #gls("lda", long: true), which is equivalent to cut $E_"xc"$ at the first order of a functional Taylor expansion: @alfeCrystalStructureThermodynamic

$ // Eq. 8.42 Alfè, Eq. 3.46 Della Pia
  E_"xc"^"LDA" [rho] := integral_V epsilon_"xc" [rho(va(r))] rho(va(r)) dif^3 va(r)
$ <eq:lda>

where $epsilon_"xc" (rho(va(r)))$ represents the exchange-correlation energy per particle of a homogeneous electron gas of density $rho(va(r))$.
In @eq:lda, the functional dependence has been substituted by a function dependence on the density, because $v_"xc" (va(r))$ has been assumed to depend only on the value of the density at the point $va(r)$:

$
  v_"xc"^"LDA" (va(r)) = dv(, rho) [epsilon_"xc" (rho(va(r))) rho(va(r))].
$

Although @lda is based on the assumption that the real density of the interacting electron system is a slowly varying function in space, it performs satisfactory well for many materials.

One limitation of @lda is that it is blind to derivatives of the charge density.
Improvements have been proposed by developing functionals that depend not just on $rho$ but also on $grad rho$.
The functionals based on this idea are known as #glspl("gga", long: true).
They are often an improvement over the @lda, but there are cases where the @lda still performs better than @gga functionals.

A more serious limitation of the @lda, or indeed of any approximation based on a sum of _local_ contributions (including @gga functionals), i.e. terms that only depend on the value of the density at the point where they are calculated, is related to their _short-sightness_.
They cannot deal with long range interactions, such as #gls("vdw", long: true) or any other electrostatic interaction that is not already encoded in the Coulomb term.

In the particular case of London dispersive interactions, which arise from the coupling of dynamically induced dipoles, appearing when a spontaneous charge fluctuation in one region of space induces the appearance of a charge fluctuation in a different region of space, information coming only from the static value of the density is not sufficient, but one also needs to relate to how the charge density changes with time in response to external pertubations.

== #md-title <sec:md>

#gls("md", long: true) is a method that allows us to sample the phase space of a isolated system of $N$ interacting classical particles, obeying the Newton's equations of motion,

$ // Eq. 7.89 Alfè
  va(f_i) = M dot.double(va(r))_i
  = - pdv(U({va(r)}), va(r_i)),
$ <eq:verlet-newton>

under certain chosen conditions, and estimate an observable expectation value (which should be an ensemble average) through the means of a time average.

Integrating @eq:verlet-newton over time, one obtains a trajectory, ${va(r)(t)}$, represented by a collection of positions of all positions of all particles in the system as a function of time.
Since the forces are conservative, the total energy of the system (kinetic plus potential) is constant and the trajectory samples the microcanonical ensemble.
If, after a sufficiently long time, the trajectory passes within an arbitrary distance from _any_ state of the ensemble and with equal probability, then the system is said to be _ergodic_ and time averages are equivalent to ensemble averages.
Ergodicity is not always satisfied, and so care must be exercised when assessing the suitability of this method for sampling the phase space.
In particular, a harmonic system is a clear example of a system that is not ergodic, as energy cannot be transferred between different normal modes.
In the following section, we discuss how to generate a trajectory in a numerical simulation, through discretization of Newton's equations of motion, @eq:verlet-newton.

=== The Verlet algorithm
The Verlet algorithm is a technique to generate the trajectory of interacting particles obeying the Newton's equations of motion, @eq:verlet-newton, through their discretization. @alfeNotesStatisticalComputational2023
Let us consider the Taylor expansion of the position of particle $i$ at time $t$, $va(r_i) (t)$, computed with forward and backward differences:

$
  arrow(r)_i (t + delta t)
  = arrow(r)_i (t) + dot(arrow(r))_i (
    t
  ) delta t + 1 / 2 dot.double(arrow(r))_i (t) (
    delta t
  )^2 + 1 / (3!) dot.triple(arrow(r))_i (t) (delta t)^3 + cal(O)((delta t)^4),
$

$
  arrow(r)_i (t - delta t)
  = arrow(r)_i (t) - dot(arrow(r))_i (
    t
  ) delta t + 1 / 2 dot.double(arrow(r))_i (t) (
    delta t
  )^2 - 1 / (3!) dot.triple(arrow(r))_i (t) (delta t)^3 + cal(O)((delta t)^4),
$

where $delta t$ is a small time interval.
Summing the two equations side by side, we obtain:
$
  arrow(r)_i (t + delta t) + arrow(r)_i (t - delta t)
  = 2 arrow(r)_i + dot.double(arrow(r))_i (t) (delta t)^2 + cal(O)((delta t)^4)
$
Consider the expression of $dot.double(arrow(r))_i$ in terms of $arrow(f)_i$ from @eq:verlet-newton, $dot.double(arrow(r))_i (t) = 1/M arrow(f)_i (t)$; substituting, we obtain:
$
  arrow(r)_i (t + delta t)
  = 2 arrow(r)_i (t) - arrow(r)_i (t - delta t) + 1 / M arrow(f)_i (t) (
    delta t
  )^2 + cal(O)((delta t)^4)
$ <eq:verlet-algorithm>
@eq:verlet-algorithm is known ad the Verlet algorithm.
@verletComputerExperimentsClassical1967[Eq. (4)]
We can re-express the equation in terms of the velocities:
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
  + cal(O)((delta t)^4)
$
$ // Eq. 7.95 Alfè
  arrow(r)_i (t + delta t)
  = arrow(r)_i (t)
  + arrow(v)_i (t) delta t
  + 1 / (2M) arrow(f)_i (t) (delta t)^2
  + cal(O)((delta t)^4),
$
which gives access to the positions at time $t + delta t$ with just the knowledge of positions, velocities, and forces at time $t$.
This expression is particularly useful at the beginning of the simulation, where only the initial positions are available, and shows that to begin a simulation we also need to provide the initial velocities.

Because of the equipartition theorem, the temperature of the system can be obtained from the ensemble average of the kinetic energy, given by:
$
  expval(E_k) = (3N) / 2 k_B T,
$
where $expval(E_k)$ is the time average of the instantaneous kinetic energy, $E_k (t)$, given by:
$
  E_k (t) = sum_i 1 / 2 M v_i^2 (t).
$

=== Ensemble averages

If the system is ergodic, ensemble averages of any physical quantity $A$ can be computed as time averages over a molecular dynamics simulation, and can be approximated as:
$
  expval(A) tilde.eq 1 / M sum_(n=1)^M A(n delta t),
$
where $A(n delta t)$ is the value of $A$ evaluated with the particles at positions ${va(r)(n delta t)}$.
The root mean square of the fluctuations of $A$ is:
$
  sigma(A) = [expval(A^2) - expval(A)^2]^(1 / 2),
$
with
$
  expval(A^2) = 1 / M sum_(n=1)^M [A(n delta t)]^2.
$

If all $M$ samples of $A$ were statistically independent from each other, then the standard deviation of $expval(A)$ would be obtained as:
$
  sigma(expval(A)) = (sigma(A)) / sqrt(M).
$
Howerver, on the trajectory generated by the molecular dynamics procedure, most of the evaluations of $A(n delta t)$ will be similar to each other, because the necessity of making small time step increments to generate an accurate trajectory means that the configurations will be close to each other for some time.
So, in computing the averages over $M$ time steps, we do not have $M$ independent samples, but only a fraction, $M/M_c$, where $M_c$ is a correlation term which measures the number of steps that we need to wait to obtain a statistically independent sample of $A$.
It follows that the statistical error on the average value of $A$ falls off as:
$
  sigma(expval(A)) = (sigma(A)) / sqrt(M/M_c).
$
$M_c$ is not known in advance, but it is fixed, and depends only on the system and on our choice of $delta t$; so, the statistical error on the average value of $A$ can be reduced as much as wanted by increasing the number of samples $M$.

=== Reblocking
// §7.3 Alfè, §4.2.1 Della Pia

The _reblocking_ procedure is a common approach to obtain $M_c$.
Suppose we split our simulation into $N$ blocks of length $M/N$ and consider the averages:
$
  expval(A)_i^N = 1 / (M / N) sum_(n = M / N i + 1)^(M / N (i+1)) A(
    n delta t
  ), quad i = 0,1,dots,N-1.
$
The average of $A$ over the whole simulation is obviously unaffected by this reblocking procedure:
$
  expval(A) = 1 / M sum_(n=1)^M A(
    n delta t
  ) = 1 / N sum_(i=0)^(N-1) expval(A)_i^N.
$
Now consider the root mean square fluctuations of the averages $expval(A)_i^N$:
$
  sigma(expval(A)^N) = [
    1 / N sum_(i=0)^(N-1) (expval(A)_i^N)^2 - expval(A)^2
  ]^(1 / 2).
$
If the evaluations of $A$ are all statistically independent from each other, then the standard deviation on the average can be obtained as:
$
  sigma_N (expval(A)) = (sigma(expval(A)^N)) / sqrt(N),
$ <eq:reblocking-sigma-N>
independently on the value of $N$, as long as $N$ is large enough to have a sufficient number of samples.
(In the extreme case of just one sample, i.e. $N=1$, it would be impossible to compute any root mean square fluctuation.)
By contrast, if the evaluations of $A$ are correlated we have:
$
  sigma_N (expval(A)) < sigma(expval(A)).
$ <eq:reblocking-inequality>
If we progressively increase the length $M/N$ of each block, at some point the averages $expval(A)_i^N$ become statistically independent.
When this happens, the inequality in @eq:reblocking-inequality becomes an equality.
Therefore, by computing $sigma_N (expval(A))$ using @eq:reblocking-sigma-N for a set of block lenghts we have a procedure to estimate $M_c$ and $sigma(expval(A))$.
Starting with $N=M$, i.e. block lenght equal to one, and working our way up reducing $N$, we find that $sigma_N (expval(A))$ increases, and at some point it reaches a plateau.
The value of $sigma_N (expval(A))$ on the plateau provides an estimate of $sigma(expval(A))$, and the block lenght corresponding to the onset of the plateau provides an estimate for the correlation length, $M_c$.

For too small values of $N$, we will observe strong oscillations in the value of $sigma_N (expval(A))$, because we are using too few evaluations of $expval(A)_i^N$ to obtain a good estimation of the error, i.e. there is no sufficient statistics.

=== Thermostats
For a sufficiently large system, averages computed in the microcanonical ensemble, with fixed $(N,V,E)$, are not much different from those computed in the canonical ensemble, with fixed $(N,V,T)$.
However, it may be desirable to be able to generate @md trajectories that span the canonical ensemble, for example because one may want to control the temperature exactly.
Moreover, some system may not be ergodic, with the harmonic system being an egregious example.
This means that given an initial set of positions and momenta, $({arrow(r)^0}, {arrow(p)^0})$, solving the Newton's equation of motion generates a trajectory that does not visit every neighbourhood of configurational space.
In such a situation time averages are biased, and do not provide good approximations for ensemble averages.

One way to overcome this problem is to couple the simulated system with an external heat bath, provided that all degrees of freedom are interacting with the bath.
Several techniques have been developed, but not all of them are capable of overcoming the ergodicity problem.
One that does is the thermostat developed by Andersen, @andersenMolecularDynamicsSimulations1980 which is based on the concept of stochastic collisions.
We know that in a perfect gas at some temperature $T$ the velocities are distributed according to the Maxwell distribution
$
  f(v) dif v
  = (m / (2 pi k_B T))^(3 / 2) 4 pi v^2 exp(-(m v^2)/(2k_B T)) dif v.
$
Therefore, one way to generate the canonical ensemble is to perform a @md simulation by repeatedly drawing velocities from a Maxwell distribution.
This periodic velocity re-initialization procedure also redistributes energy between different modes, and so it is an effective way to overcome the ergodicity problem.
It can be shown that the frequency of these velocity re-initializations does not affect the ability to sample the canonical ensemble; however, drawing the velocities too often will result in the system moving very slowly from one region of configuration space to another; drawing them too seldom results in slow transfer of energy between different modes, which would only overcome the ergodicity problem slowly.
Finding the appropriate time interval between velocities randomizations is then a matter of finding the right compromise to maximize efficiency.

// https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.langevin
In *Langevin dynamics*#footnote[https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.langevin], a (small) friction term and a fluctuating force are added to Newton's second law, @eq:verlet-newton, which is then integrated numerically.
The temperature of the heat bath and magnitude of the friction is specified by the experimenter, the amplitude of the fluctuating force is then calculated to give that temperature.
This procedure has some physical justification: in a real metal the atoms are (weakly) coupled to the electron gas, and the electron gas therefore acts like a heat bath for the atoms.
If heat is produced locally, the atoms locally get a temperature that is higher than the temperature of the electrons, heat is transferred to the electrons and then rapidly transported away by them.
A Langevin equation is probably a reasonable model for this process.

A disadvantage of using Langevin dynamics is that if significant heat is produced in the simulation, then the temperature will stabilize at a value higher than the specified temperature of the heat bath, since a temperature difference between the system and the heat bath is necessary to get a finite heat flow.
Another disadvantage is that the fluctuating force is stochastic in nature, so repeating the simulation will not give exactly the same trajectory, if not using exactly the same starting configuration and random number generator.

In *Andersen dynamics*#footnote[https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.andersen], constant temperature is imposed by stochastic collisions with a heat bath.
With a (small) probability the collisions act occasionally on velocity components of randomly selected particles.
Upon a collision the new velocity is drawn from the Maxwell-Boltzmann distribution at the corresponding temperature.
The system is then integrated numerically at constant energy according to the Newtonian laws of motion.
The collision probability is defined as the average number of collisions per atom and timestep.
The algorithm generates a canonical distribution. @daanfrenkelUnderstandingMolecularSimulation2002[§6.1.1]
However, due to the random decorrelation of velocities, the dynamics are unphysical and cannot represent dynamical properties like e.g. diffusion or viscosity.
Another disadvantage is that the collisions are stochastic in nature, so repeating the simulation will not give exactly the same trajectory, if not using exactly the same starting configuration and random number generator.
Typical values for the collision probability are in the order of $10^(-4) ÷ 10^(-1)$.

In *Nosé-Hoover dynamics*, an extra term is added to the Hamiltonian representing the coupling to the heat bath.
From a pragmatic point of view one can regard Nosé-Hoover dynamics as adding a friction term to Newton's second law, but dynamically changing the friction coefficient to move the system towards the desired temperature.
Typically the "friction coefficient" will fluctuate around zero.

During simulations in the present work, the Langevin thermostat was used for constant $(N,V,T)$ @md and combined Nose-Hoover and Parrinello-Rahman#footnote[https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.npt] dynamics for the $(N, P, T)$ ensemble.
The Berendsen thermostat was considered, but later discarded#footnote[See the "Flying ice cube" effect.] in favour of the thermostats above.

=== Mean square displacement
// §7.3.3 Alfè, §4.2.3 Della Pia
In a molecular dynamics simulation, a convenient quantity that can be used to monitor the state of the system is the _mean square displacement_, defined as:
$
  m(t) := 1 / N sum_(i=1)^N |va(r_i)(t + t_0) - va(r_i)(t_0)|^2,
$ <eq:mean-square-displacement>
where $t_0$ is some initial reference time.
In a system with no diffusing behaviour, such as a solid, $m(t)$ is expected to rise first, and then reach a constant, which is related to the maximum displacement of the particles from their equilibrium positions.
By contrast, in a fluid $m(t)$ is expected to rise with time.
If the motion of the particles is random, which is a good approximation for a system in thermal equilibrium, then $m(t)$ increases linearly with time, and its slope is related to the _diffusion coefficient_.

To understand where the linear behaviour of $m(t)$ comes from, consider a random walk, for simplicity of a one-dimensional system.
Let us associate a variable $z_i$ to the $i^"th"$ step in the walk, which can be either $+1$ or $-1$, depending if the step is taken by going to the right or to the left.
Let us also define the variable $s_N := sum_(i=1)^N z_i$, which is the length of the walk after $N$ steps.
The average value of $s_N$ is zero, as:
$
  expval(s_N) = expval(sum_(i=1)^N z_i) = sum_(i=1)^N expval(z_i) = 0.
$
The average of $s_N^2$, however, is not zero:
$
  expval(s_N^2)
  = expval(sum_(i=1)^N z_i sum_(j=1)^N z_j)
  = sum_(i,j=1)^N expval(z_i z_j)
  = sum_(i=1)^N underbrace(expval(z_i^2), 1) + sum_(i,j=1 \ i != j)^N underbrace(expval(z_i z_j), expval(z_i) expval(z_j) = 0)
  = N,
$
as $expval(z_i^2) = 1$, and $expval(z_i z_j) = expval(z_i) expval(z_j) = 0 dot 0 = 0$, since the $i^"th"$ and the $j^"th"$ steps are uncorrelated.
This shows that in a random walk of step size 1 the mean square displacement from the origin of the walk is equal to the number of steps, and so it is linearly proportional to time, if the number of steps per unit time is constant.
In systems with continuous displacements, this translates into a linear dependence on time of the mean square displacement $m(t)$ defined above.

Over a simulation of total length $T$, one clearly only has access to $m(t)$ with $0 <= t <= T$, and to improve on statistics it is useful to compute @eq:mean-square-displacement by averaging over time origins $t_0$:
$
  m(t)
  = 1 / (T - t) sum_(t_0=0)^(T-t) 1 / N sum_(i=1)^N |va(r_i)(t+t_0) - va(r_i)(
    t_0
  )|^2.
$
We see that for $t=0$ it is possible to average on the whole length of the simulation, but, as $t$ increases, the available length over which one can average is reduced to $T-t$, and so the statistical error on $m(t)$ increases with $t$.
For $t=T$ there is only one available configuration.

The mean square displacement becomes particularly useful in @md simulations at high temperature, i.e. in proximity of the melting temperature, to check if the system is in the solid or liquid phase.
This check is also important in thermodynamic integration simulations with all contributions included ($lambda = 1$), where the potential is in the complete form, $U$, and there is no explicit harmonic contribution that makes atoms oscillate around their equilibrium positions.

== #thermodynamics-crystal-stability-title <sec:thermodynamics-crystal-stability>

The main quantity to consider to assess the stability of a crystal is the lattice free energy, which is defined as the free energy per molecule gained upon assuming the crystal form with respect to the gas state.
This is a very expensive quantity to estimate because it is a finite temperature property.
For this reason, we often approximate finite temperature free energy differences with zero temperature lattice energy, $E_"latt"$, which is analogously defined as the energy per molecule gained upon assuming the crystal form with respect to the gas state.
It can be computed as:
$
  E_"latt" = E_"crys" - E_"gas",
$ <eq-zen_2018_1>
with $E_"crys"$ as the total electronic energy per molecule in the crystal state and $E_"gas"$ as the energy of the isolated molecule.

Typically, the computation of $E_"latt"$ is performed at zero temperature and considering only the electronic contribution; i.e., zero-point motion and quantum nuclear effects are neglected. @beranPredictingMolecularCrystal2016
The lattice energy is not directly assessable experimentally, but it can be indirectly obtained from experimental measures of the sublimation enthalpy, $Delta_"sub" H(T)$, at a given temperature $T$, by including a (theoretically evaluated) energy contribution, $Delta_"T&QN" (T)$, accounting for thermal and quantum nuclear effects:
$
  Delta_"sub" H(T) = -E_"latt" + Delta_"T&QN" (T).
$
The evaluation of $Delta_"T&QN" (T)$ can be challenging, especially for large molecules where anharmonic contributions are important. @reillyUnderstandingRoleVibrations2013
Since both $Delta_"sub" H(T)$ and $Delta_"T&QN" (T)$ are affected by errors, accurate theoretical evaluations of $E_"latt"$ are of help for comparison.

In order to derive $Delta_"T&QN"$, we need to start from the definition of the sublimation enthalpy, $Delta_"sub" H(T)$, that is the difference between the enthalpy of the gas, $H^g (T)$, and of the crystal solid, $H^s (T)$, both at temperature $T$.
By separating the electronic ($"el"$), translational ($"trans"$), rotational ($"rot"$) and vibrational ($"vib"$) contributions, and noticing that in the crystal there are no trans-rotational contributions, we have that
$
  Delta_"sub" H = E_"el"^g + E_"trans"^g + E_"rot"^g + E_"vib"^g + p V - (
    E_"el"^s + E_"vib"^s
  ),
$ <eq-zen_si_14>
where the superscript stands either for gas ($g$) or solid ($s$), and the temperature dependance has been dropped for the seek of brevity.
We will also consider as negligible the pressure times volume term, $p V$.
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

Dispersion interactions, sometimes called van der Waals interactions, are crucial for describing the weak, long-range interactions between electrons.
For many systems it is important to take dispersion interactions into account for reliable prediction,
as they are essential in a collection of active research fields in solid-state physics and chemistry, including molecular crystal packing, crystal structure prediction, surface adsorption and reactivity, and supramolecular chemistry.
The representation of dispersion interactions is not possible within common approximations in @dft, like local or semi-local functionals such as PBE, because dispersion arises from non-local correlation effects involving distant fragments in the crystal. @otero-de-la-rozaBenchmarkNoncovalentInteractions2012

This motivates the use of additive non-local corrections.
In @dft, this is typically done either by applying an a posteriori correction to a certain @dft functional prediction, $E_"DFT"$, or by including non-local terms in the exchange-correlation functional.
In the first of the two cases, the final energy is computed as
$ E = E_"DFT" + E_"disp", $ <eq-otero_2012_1>
where $E_"disp"$ is a function of the electron-electron distance with parameters that are fitted for the specific functional the correction is applied to. (For instance, the D3 parameters for @pbe will be different from the D3 parameters of revPBE.)
Inclusion of a dispersion correction to @dft is necessary to describe the dynamics of liquid water, the geometries and binding energies of layered solids, and stability of metal-organic frameworks, among many other examples. @batatiaFoundationModelAtomistic2023[§4]

Several alternatives for the correction $E_"disp"$ are available (D2, D3, XDM, D4, TS, MBD, ...), but in principle we do not know which one yields the most reliable prediction for a specific system.
Additive dispersion corrections typically employ a physical model for dispersion interactions with empirical parameters optimized to cut off the correction at interatomic distances where approximate @dft is reliable.
As illustrative examples of methods to take dispersion into account, the XDM and D3 corrections are detailed below, as they are employed in the reference implementations cited for this thesis.

==== Exchange-hole Dipole Moment

The @xdm model describes the dispersion energy of two neutral fragments as the electrostatic interaction of the dipoles formed by electrons and their associated exchange holes.
Here, $E_"disp"$ contains the usual $R^(-6)$ leading term as well as two additional higher order atomic-pairwise terms
$
  E_"disp" = -1 / 2 sum_(i j) sum_(n = 6,8,10) (C_(n,i j)) / (
  R^n_("vdw", i j) + R^n_(i j)
  ).
$ <eq-otero_2012_2>
The fundamental objects in this equation are the inter-atomic interaction coefficients $C_(n, i j)$ that in the @xdm model are calculated entirely from first-principles quantities using second-order perturbation theory. @otero-de-la-rozaBenchmarkNoncovalentInteractions2012[Equations 3-5]

All the objects above are parameter-free, except for the damping expression in @eq-otero_2012_2.
The interatomic van der Waals radii ($R_("vdw", i j)$) are a set of parametric values that control the distance at which the pairwise dispersion interactions are switched off.
The value of the parameters in $R_("vdw", i j)$ is obtained by fitting to a training set both in gas-phase and under periodic-boundary conditions. @otero-de-la-rozaBenchmarkNoncovalentInteractions2012[Equations 8, 9]

Because the dispersion coefficients are calculated rather than ﬁtted,
@eq-otero_2012_1 works under the assumption that the @dft functional presents a completely dispersionless behavior.
This requirement is not met by most @gga functionals,
which are sometimes too repulsive and sometimes spuriously binding,
depending on the reduced-density-gradient tail behavior of the exchange enhancement factors.

==== DFT-D3
DFT-D3 is an interatomic potential which uses tabulated values of atomic polarizabilities to describe two-body and, optionally, three-body Axilrod-Teller dispersion interactions.
As MACE-MP-0 is trained to @pbe energies, forces, and stresses, it inherits @pbe's lack of long-range dispersion interactions.
An optional, additive DFT-D3 dispersion correction can be applied to MACE-MP-0.
The same parameters used in PBE-D3(BJ), i.e., DFT-D3 with a Becke-Johnson damping function @grimmeEffectDampingFunction2011, are used in the D3 correction to MACE-MP-0. @batatiaFoundationModelAtomistic2023[§4]

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

// === The role of anharmonic contributions
// // TODO
// @rossiAnharmonicQuantumFluctuations2016

// The harmonic treatment ignores effects due to anharmonic thermal motion and cell expansion. @reillyUnderstandingRoleVibrations2013

// In the harmonic approximation, supercell phonon calculations show signiﬁcant deviations from the widely used Dulong-Petit law, as noted elsewhere. @reillyUnderstandingRoleVibrations2013

// Assessing predictions of lattice energies requires careful consideration of vibrational, many-body dispersion and exact-exchange contributions. @reillyUnderstandingRoleVibrations2013

// Exact exchange, which is rarely considered in @dft studies of molecular crystals, is shown to have a signiﬁcant contribution to lattice energies, systematically improving agreement between theory and experiment. @reillyUnderstandingRoleVibrations2013

// Hybrid functionals are often not used in the study of cohesive properties of molecular crystals,
// @al-saidiAssessmentVdWTSMethod2012, @otero-de-la-rozaBenchmarkNoncovalentInteractions2012
// largely due to their additional computational cost, which particularly in a plane-wave basis can reach more than an order of magnitude larger than the corresponding semi-local functional. @reillyUnderstandingRoleVibrations2013

// #box(
//   stroke: 2pt + red,
//   inset: 1mm,
//   [The pharmaceutical industry spends considerable resources on high-throughput crystallization experiments to screen for polymorphs, // TODO citazione a S. L. Morissette et al., High-throughput crystallization: Polymorphs, salts, co-crystals and solvates of pharmaceutical solids. Adv. Drug Deliv. Rev. 56, 275–300 (2004).
//     into which the target structure may decay.
//     However, crystallization experiments do not probe thermodynamic stability, and conclusive studies of the impact of temperature changes after crystallization on the stability of polymorphs (i.e. their monotropic or enantiotropic nature) // TODO cite E. H. Lee, A practical guide to pharmaceutical polymorph screening & selection. Asian J. Pharm. Sci. 9, 163–175 (2014)
//     are often prevented by limited sample quantities.
//     Hence, there is the appeal of theoretical @csp // TODO cite 14 S. L. Price, Predicting crystal structures of organic compounds. Chem. Soc. Rev. 43,2098–2111 (2014).
//     based on the thermodynamic stability, which promises to complement crystallization experiments // TODO cite 15 J. Nyman, S. M. Reutzel-Edens, Crystal structure prediction is changing from basic science to applied technology. Faraday Discuss. 211, 459–476 (2018).
//     by exhaustively searching for competing polymorphs.
//     @kapilCompleteDescriptionThermodynamic2022
//   ],
// )

// #box(
//   stroke: blue + 2pt,
//   inset: 1mm,
//   [
//     The size and sign of nuclear quantum effects, anharmonicity, and cell expansion and flexibility, depend entirely on the compound and the polymorphs at hand, highlighting that rigorous @qti is indispensable for predicting phase stabilities and that molecular crystals are typically stabilized by a nontrivial interplay of different physical effects, whose individual importance is belied by the subtle resultant free energy differences. @kapilCompleteDescriptionThermodynamic2022
//   ],
// )

== #machine-learning-title <sec:machine-learning>
#text(blue)[
  Machine Learning can be described as _the application and science of algorithms that make sense of data_. @raschkaMachineLearningPyTorch2022[§1]
  There are three types of machine learning: supervised learning, unsupervised learning, and reinforcement learning.
  It is our interest to study the *supervised learning* type.
  The main goal in supervised learning is to learn a model from labeled training data that allows us to make predictions about unseen or future data.
  Here, the term "supervised" refers to a set of training examples (data inputs) where the desired output signals (labels) are already known.
  Supervised learning is then the process of modeling the relationship between the data inputs and the labels.
  Thus, we can also think of supervised learning as "label learning".
  A supervised learning task with discrete class labels is also called a *classification task*. Another subcategory of supervised learning is *regression*, where the outcome signal is a continuous value.
]
In the detailed description of MACE in @sec:mace, we will see that the data inputs are atom positions, atomic number; while the label will be the energy, a continuous value learned through regression.
#text(blue)[
  In regression analysis, we are given a number of predictor (*explanatory*) variables and a continuous response variable (*outcome*), and we try to find a relationship between those variables that allows us to predict an outcome.
  In the field of machine learning, the predictor variables are commonly called "features", and the response variables are usually referred to as "target variables".
]

#text(blue)[
  Artificial neurons represent the building blocks of the multilayer artificial #glspl("nn").
  The basic concept behind artificial #glspl("nn") was built upon hypotheses and models of how the human brain works to solve complex problem tasks.
  NNs are more popular today than ever thanks to the many breakthroughs that have been made in the previous decade, which resulted in what we now call deep learning algorithms and architectures—NNs that are composed of many layers.

  === Single layer neural network
  Before we dig deeper into a particular multilayer @nn architecture, let’s briefly reiterate some of the concepts of single-layer NNs, namely, the @adaline algorithm.

  #figure(
    image("thesis/imgs/raschkaMachineLearningPyTorch2022_11_01.png"),
    caption: [
      The Adaline algorithm. Taken from @raschkaMachineLearningPyTorch2022[Fig. 11.1].
    ],
  )

  The Adaline algorithm performs binary classification, and uses the gradient descent optimization algorithm to learn the weight coefficients of the model.
  In every epoch (pass over the training dataset), we update the weight vector $va(w)$ and bias unit $b$ using the following update rule:
  $
    va(w) :=
    va(w) + Delta va(w), quad
    b := b + Delta b,
  $
  where $Delta w_j = - eta pdv(L, w_j)$ for each weight $w_j$ in the weight vector $va(w)$ and $Delta b = - eta pdv(L, b)$ for the bias unit.
  In other words, we computed the gradient based on the whole training dataset and updated the weights of the model by taking a step in the opposite direction of the loss gradient $grad L (va(w))$.
  In order to find the optimal weights of the model, we optimize an objective function that we define as the mean of squared errors (MSE) loss function $L(va(w))$.
  Furthermore, we multiply the gradient by a factor, the learning rate $eta$, which we have to choose carefully to balance the speed of learning against the risk of overshooting the global minimum of the loss function.
  In gradient descent optimization, we update all weights simultaneously after each epoch, and we define the partial derivative for each weight $w_j$ in the weight vector, $va(w)$, as follows:
  $
    pdv(L, w_j) =
    pdv(, w_j) 1 / n sum_i (y^((i)) - a^((i)))^2 = - 2 / n sum_i (
      y^((i)) - a^((i))
    ) x_j^((i)).
  $
  Here, $y^((i))$ is the target class label of a particular sample $x^((i))$, and $a^((i))$ is the activation of the neuron, which is a linear function in the special case of Adaline.
  Furthermore, one can define the activation function $sigma(dot)$ as follows:
  $
    sigma(dot) = z = a.
  $
  Here, the net input, $z$, is a linear combination of the weights that are connecting the input layer to the output layer:
  $
    z =
    sum_j w_j x_j + b = va(w)^TT va(x) + b.
  $
  While the activation $sigma(dot)$ is used to compute the gradient update, a threshold function is implemented to squash the continuous-valued output into binary class labels for prediction:
  $
    hat(y) = cases(
      1 & "if" z >= 0,
      0 & "otherwise"
    )
  $
  A frequent technique used to accelerate the model training is the so-called *stochastic gradient descent (SGD)* optimization.
  SGD approximates the loss from a single training example (online learning) or a small subset of training examples (mini-batch learning).
  Apart from faster learning---due to the more frequent weight updates compared to gradient descent---its noisy nature is also regarded as beneficial when training multilayer NNs with nonlinear activation functions, which do not have a convex loss function.
  Here, the added noise can help to escape local loss minima.

  === The multilayer neural network architecture
  Here we will explain how to connect multiple single neurons to a multilayer feedforward @nn; this special type of _fully connected_ network is also called Multi-Layer Perceptron (MLP).

  #figure(
    image("thesis/imgs/raschkaMachineLearningPyTorch2022_11_02.png"),
    caption: [A two-layer MLP. Figure from @raschkaMachineLearningPyTorch2022.],
  ) <fig:multi-layer-perceptron>

  Next to the data input, the MLP depicted in @fig:multi-layer-perceptron has one hidden layer and one output layer.
  The units in the hidden layer are fully connected to the input features, and the output layer is fully connected to the hidden layer.
  If such a network has more than one hidden layer, we alco call it a *deep NN*.
  The number of layers and units in a NN can be tought of as additional hyper-parameters that we want to optimize for a given problem.

  We denote the $i$th activation unit in the $l$th layer as $a_i^((l))$.
  We use the "in" superscript for the input features, the "h" superscript for the hidden layer, and the "out" superscript for the output layer.
  The $va(b)$'s denote bias units; those are vectors with the number of elements being equal to the number of nodes in the layer they correspond to.

  Each node in layer $l$ is connected to all nodes in layer $l + 1$ via a weight coefficient.
  The connection between the $k$th unit in layer $l$ to the $j$th unit in layer $l + 1$ will be written as $w_(j,k)^((l))$.
  We denote the weight matrix that connects the input to the hidden layer as $W^(("h"))$, and we write the matrix that connects the hidden layer to the output layer as $W^(("out"))$.

  The usage of multiple units in the output layer allow for native multi-class classification, via a generalization of the one-versus-all (OvA) technique.
  This is similar to the one-hot representation of categorical variables.

  To calculate the output of a MLP model, we employ the process of *forward propagation*.
  The MLP learning procedure can be summarized in three steps:

  + Starting at the input layer, forward propagate the patterns of the training data through the network to generate an output.
  + Based on the network's output, calculate the loss that we want to minimize using a loss functions of choice.
  + Backpropagate the loss, find its derivative with respect to each weight and bias unit in the network, and update the model.

  Finally, after we repeat these three steps for multiple epochs and learn the weight and bias parameters of the MLP, we use forward propagation to calculate the network output.

  Since each unit in the hidden layer is connected to all units in the input layer, we first calculate the activation unit of the hidden layer $a_1^((h))$ as follows:
  $
    z_1^((h)) &= x_1^(("in")) w_(1,1)^((h)) + x_2^(("in")) w_(1,2)^((
      h
    )) + dots + x_m^(("in"))w_(1,m)^((h)), \
    a_1^((h)) &= sigma(z_1^((h))).
  $
  Here, $z_1^((h))$ is the net input and $sigma(dot)$ is the activation function, which has to be differentiable to learn the weights that connect the neurons using a gradient-based approach.
  To be able to solve complex problems, the MLP model needs nonlinear activation functions, for example, the sigmoid (logistic) activation function:
  $
    sigma(z) = 1 / (1 + e^(-z))
  $

  #figure(
    image("thesis/imgs/raschkaMachineLearningPyTorch2022_11_03.png"),
    caption: [
      The sigmoid activation function.
      Figure from @raschkaMachineLearningPyTorch2022.
    ],
  )

  The MLP is a typical example of a feedforward artificial NN.
  The term *feedforward* refers to the fact that each layer serves as the input to the next layer without loops, in contrast to recurrent NNs.

  The activation is usually written in a more compact, vectorized form, using the concepts of linear algebra:
  $
    va(z)^((h)) &= va(x)^(("in")) W^((h)TT) + va(b)^((h)), \
    va(a)^((h)) &= sigma(va(z)^((h))).
  $
  Here, $va(z)^((h))$ is a $1 times m$ dimensional feature vector;
  $W^((h))$ is a $d times m$ dimensional weight matrix, where $d$ is the number of units in the hidden layer;
  the bias vector $va(b)^((h))$ consists of $d$ bias units (one bias unit per node).

  Generalizing the computation to all $n$ samples in the training dataset:
  $
    Z^((h)) = X^(("in")) W^((h)TT) + va(b)^((h))
  $
  Here, $X^(("in"))$ is a $n times m$ matrix, and the matrix multiplication will result in a $n times d$ dimensional net input matrix, $Z^((h))$.
  Finally, we apply the activation function $sigma(dot)$ to each value in the net input matrix to get the $n times d$ dimensional activation matrix in the next layer:
  $
    A^((h)) = sigma(Z^((h)))
  $
  Similarly for the output layer we get:
  $
    Z^(("out")) = A^((h)) W^(("out") TT) + va(b)^(("out")), quad
    A^(("out")) = sigma(Z^(("out"))).
  $

  === Training an artificial neural network

  We use an MSE loss (as in Adaline) to train the multilayer NN as it makes the derivation of the gradients a bit easier to follow.
  If we predict the class label of an input data with class label 2, using this MLP, the activation of the third layer and the target may look like this:
  $
    a^(("out")) =
    vec(0.1, 0.9, dots.v, 0.3), quad
    y = vec(0, 1, dots.v, 0).
  $

  Thus, our MSE loss either has to sum or average over the $t$ activation units in our network in addition to averaging over the $n$ examples in the dataset or mini-batch:
  $
    L(W, va(b)) =
    1 / n sum_1^n 1 / t sum_(j=1)^t (y_j^((i)) - a_j^(("out")(i)))^2
  $

  The goal is to minimize the loss function $L(W)$.
  We need to calculate the partial derivative of the parameters $W$ with respect to each weight for every layer in the network:
  $
    pdv(L, w_(j,l)^((i)))
  $

  Note that $W$ consists of multiple matrices.
  In an MLP with one hidden layer, we have the weight matrix, $W^((h))$, which connects the input to the hidden layer, and $W^(("out"))$, which connects the hidden layer to the output layer.

  Backpropagation is a very efficient and one of the most widely used algorithms for training artificial NNs.
  In essence, we can think of backpropagation as a very computationally efficient approach to compute the partial derivatives of a complex, non-convex loss function in multilayer NNs.
  The goal is to use those derivatives to learn the weight coefficients for parameterizing such a multilayer artificial NN.
  The error surface of an NN loss function is not convex or smooth with respect to the parameters.
  There are many bumps in this high-dimensional loss surface (local minima) that we have to overcome in order to find the global minimum of the loss function.
  Given the function $F(x) = f(g(h(u(v(x)))))$, we can use the chain rule to compute the derivative:
  $
    dv(F,x) =
    dv(,x) f(g(h(u(v(x))))) =
    dv(f,g) dv(g,h) dv(h,u) dv(u,v) dv(v,x).
  $
  In the context of computer algebra, a set of techniques, known as *automatic differentiation*, has been developed to solve such problems very efficiently.

  Automatic differentiation comes with two modes, the forward and reverse modes; backpropagation is simply a special case of reverse-mode automatic differentiation.
  The key point is that applying the chain rule in forward mode could be quite expensive since we would have to multiply large matrices for each layer (Jacobians) that we would eventually multiply by a vector to obtain the output.

  The trick of reverse mode is that we traverse the chain rule from right to left.
  We multiply a matrix by a vector, which yields another vector that is multiplied by the next matrix, and so on.
  Matrix-vector multiplication is computationally much cheaper than matrix-matrix multiplication, which is why backpropagation is one of the most popular algorithms used in NN training.

  #figure(
    image("thesis/imgs/raschkaMachineLearningPyTorch2022_11_13.png"),
    caption: [
      Computing the partial derivatives of the loss with respect to the first hiddel layer weight.
      Averaging over the mini-batch is omitted.
      Figure from @raschkaMachineLearningPyTorch2022.
    ],
  )
]

== #gnn-title <sec:gnn>
#text(blue)[Neural networks have been adapted to leverage the structure and properties of graphs.]

#figure(
  image("thesis/imgs/distill.pub.gnn-intro.fig1.png"),
  caption: [
    Taken from @sanchez-lengelingGentleIntroductionGraph2021.
    Diagram representing of how a node accumulates information from nodes around it through the layers of the network.
    // Hover over a node in the diagram below to see how it accumulates information from nodes around it through the layers of the network.
  ],
)

#text(blue)[A set of objects, and the connections between them, are naturally expressed as a graph.]
#text(blue)[#glspl("gnn") see practical applications in physics simulations @sanchez-gonzalezLearningSimulateComplex2020] and are the foundation of the main calculator employed in this work, MACE, as detailed in @sec:mace.

#text(blue)[A graph represents the relations (_edges_) between a collection of entities (_nodes_).]

#large_figure(
  grid(
    columns: 3,
    gutter: 10pt,
    figure(
      image("thesis/imgs/distill.pub.gnn-intro.what-is-a-graph-V.png"),
      caption: [*V* -- Vertex (or node) attributes e.g., node identity, number of neighbors.],
      numbering: none,
      outlined: false,
    ),
    figure(
      image("thesis/imgs/distill.pub.gnn-intro.what-is-a-graph-E.png"),
      caption: [*E* -- Edge (or link) attributes and directions e.g., edge identity, edge weight.],
      numbering: none,
      outlined: false,
    ),
    figure(
      image("thesis/imgs/distill.pub.gnn-intro.what-is-a-graph-U.png"),
      caption: [*U* -- Global (or master node) attributes e.g., number of nodes, longest path.],
      numbering: none,
      outlined: false,
    ),
  ),
  caption: [
    Three types of attributes we might find in a graph.
    Figures from @sanchez-lengelingGentleIntroductionGraph2021.
  ],
)

#text(blue)[
  To further describe each node, edge or the entire graph, we can store information in each of these pieces of the graph.
  We can additionally specialize graphs by associating directionality to edges (_directed_, _undirected_).

  Graphs are very flexible data structures], and for this reason they were used by MACE to embed atomistic properties of physical systems.
#text(blue)[
  It’s a very convenient and common abstraction to describe this 3D object as a graph, e.g. where nodes are atoms and edges are covalent bonds.
]
#text(blue)[
  A way of visualizing the connectivity of a graph is through its adjacency matrix.
  One labels the nodes, in this case each of 14 non-H atoms in a caffeine molecule, and fill a matrix of $n_"nodes" times n_"nodes"$ with an entry if two nodes share an edge.
]

#figure(
  image("thesis/imgs/distill.pub.gnn-intro-caffeine-adiacency-graph.png"),
  caption: [
    Figure from @sanchez-lengelingGentleIntroductionGraph2021.
    (Left) 3D representation of the Caffeine molecule. (Center) Adjacency matrix of the bonds in the molecule. (Right) Graph representation of the molecule.
  ],
)

*Definition*:
#text(blue)[
  A GNN is an optimizable transformation on all attributes of the graph (nodes, edges, global-context) that preserves graph symmetries (permutation invariances).
]

#text(blue)[
  In the following, we will describe the @mpnn framework proposed by @gilmerNeuralMessagePassing2017 using the Graph Nets architecture schematics introduced by @battagliaRelationalInductiveBiases2018.
  #glspl("gnn") adopt a “graph-in, graph-out” architecture meaning that these model types accept a graph as input, with information loaded into its nodes, edges and global-context, and progressively transform these embeddings, without changing the connectivity of the input graph.

  The GNN uses a differentiable model of choice (e.g. a multilayer perceptron (MLP)) on each component of a graph; this is a GNN layer.
  For each node vector, one applies the model and gets back a learned node-vector.
  One does the same for each edge, learning a per-edge embedding, and also for the global-context vector, learning a single embedding for the entire graph.

  #figure(
    image("thesis/imgs/distill.pub.gnn-intro.single-layer.png"),
    caption: [
      A single layer of a simple GNN. A graph is the input, and each component (V, E, U) gets updated by a MLP to produce a new graph.
      Each function subscript indicates a separate function for a different graph attribute at the n-th layer of a GNN model.
    ],
  )

  As is common with neural network modules or layers, one can stack these GNN layers together.

  Because a GNN does not update the connectivity of the input graph, one can describe the output graph of a GNN with the same adjacency list and the same number of feature vectors as the input graph.
  But, the output graph has updated embeddings, since the GNN has updated each of the node, edge and global-context representations.

  === GNN Predictions by Pooling Information
  How does a GNN make predictions in any of its tasks described above?
  Prediction tasks can belong to binary classification, multi-class classification or regression cases.
  In the following example, the binary classification will be considerer for brevity, but this framework extends to the other cases.
  If the task is to make predictions on nodes, and the graphs already contains node information, the approach is straightforward---for each node embedding, apply a linear classifier.

  #figure(image("thesis/imgs/distill.pub.gnn-intro.linear-classifier.png"))

  However, it is not always so simple.
  For instance, one might have information in the graph stored in edges, and no information in nodes,but still need to make predictions on nodes.
  One needs a way to collect information from edges and give them to nodes for prediction.
  One can do this by _pooling_.
  Pooling proceeds in two steps:

  + For each item to be pooled, _gather_ each of their embeddings and concatenate them into a matrix.
  + The gathered embeddings are then _aggregated_, usually via a sum operation.

  @sanchez-lengelingGentleIntroductionGraph2021 represents the _pooling_ operation by the letter $rho$, and denotes that we are gathering information from edges to nodes as $p_(E_n arrow V_n)$.

  #figure(
    image("thesis/imgs/distill.pub.gnn-intro.aggregate-information-from-adjacent-edges.png"),
    caption: [The edges connected to the black node are gathered and aggregated to produce an embedding for that target node. Figure from @sanchez-lengelingGentleIntroductionGraph2021.],
  )

  So if we only have edge-level features, and are trying to predict node information, we can use pooling to route (or pass) information to where it needs to go, the model looks like this:
  #figure(
    image("thesis/imgs/prediction_edges_nodes.e6796b8e.png"),
    caption: [Figure from @sanchez-lengelingGentleIntroductionGraph2021.],
  )
  The same reasoning goes for the prediction of edge-level information and global properties, gathering available node and/or edge information together and aggregating them to get the desired predictions.
  This technique is similar to _Global Average Pooling_ layers in #glspl("cnn").
  #figure(
    image("thesis/imgs/prediction_nodes_edges_global.7a535eb8.png"),
    caption: [
      This is a common scenario for predictiong molecular properties.
      Figure from @sanchez-lengelingGentleIntroductionGraph2021.
    ],
  )
  The classification model $cal(C)$ can easily be replaced with any differentiable model, or adapted to multi-class classification using a generalized linear model.
  #figure(
    image("thesis/imgs/Overall.e3af58ab.png"),
    caption: [
      An end-to-end prediction task with a GNN model.
      Figure from @sanchez-lengelingGentleIntroductionGraph2021.
    ],
  )

  === Passing messages between parts of the graph
  One could make more sophisticated predictions by using pooling within the @gnn layer, in order to make the learned embeddings aware of graph connectivity.
  We can do this using _message passing_ @gilmerNeuralMessagePassing2017, where neighbouring nodes or edges exchange information and influence each other's update embeddings.

  Message passing works in three steps:

  + For each node in the graph, _gather_ all the neighbouring node embeddings (or messages).
    // which is the $g$ function described above.
    // <Non vedo nessuna funzione g nell'articolo>
  + Aggregate all messages via an aggregate function (like sum).
  + All pooled messages are passed through an _update function_, usually a learned neural network.

  Just as pooling can be applied to either nodes or edges, message passing can occur between either nodes or edges.

  These steps are key for leveraging the connectivity of graphs.
  One can build more elaborate variants of message passing in @gnn layers that yield @gnn models of increasing expressiveness and power.

  #figure(
    image("thesis/imgs/distill.pub.gnn-intro.message-passing-update.png"),
    caption: [Pooling, update and storage of the adjacent embedding for the highlighted node. Figure from @sanchez-lengelingGentleIntroductionGraph2021.],
  )

  This sequence of operations, when applied once, is the simplest type of message-passing @gnn layer.
  This is reminiscent of standard convolution: in essence, message passing and convolution are operations to aggregate and process the information of an element's neighbours in order to update the element's value.
  In graphs, the element is a node, and in images, the element is a pixel.
  However, the number of neighbouring nodes in a graph can be variable, unlike in an image where each pixel has a set number of neighbouring elements.

  By stacking messge passing @gnn layers together, a node can eventually incorporate information from across the entire graph: after three layers, a node has information about the nodes three steps away from it.

  The updated architecture diagram to include this new source of information for nodes is the following:
  #figure(
    image("thesis/imgs/arch_gcn.40871750.png"),
    caption: [
      Schematic for a GCN architecture, which updates node representations of a graph by pooling neighbouring nodes at a distance of one degree.
      Figure from @sanchez-lengelingGentleIntroductionGraph2021.
    ],
  )
]

The notions exposed above are a sufficient introduction to understand the basic functioning of the MACE calculator, described in detail in @sec:mace, that was used for the work in this thesis, the results of which are available in the next chapters, @sec:results-1 and @sec:results-2.

= Results I: model assessment <sec:results-1>

In this chapter we test the MACE calculators on small, known systems,
about which the literature is abundant and a direct comparison of results is possible.
We analyze the geometrical and vibrational features of the water molecule and the water dimer after optimization with each calculator.
The binding energy of the dimer is calculated and compared with reference @dft methods.
In particular, we test three versions of the foundation model MACE-MP-0, namely _small_, _medium_ and _large_.
The three models differ in the number of layers used in the MACE architecture, hence the total number of parameter built into the model.
Moreover, a fine-tuned @mlp for ice polymorphs @kaurDataefficientFinetuningFoundational2024, MACE-ICE13-1, is also put to the test and compared with the MACE foundation models.

== The water molecule <sec-molecule>

The first task is the optimization of the geometry and the calculation of vibrational properties of the water molecule.
We studied the relaxation of the geometry of a single water molecule, also referred to as the monomer, and compared the results with physical values of the gas phase molecule. The optimization is tackled with two concurring methods:
static local *minimization of the potential energy* and *analysis of vibrational properties* to assess the dynamical stability.

#figure(
  image("simulazioni/02_water/01_molecule/MACE-ICE13-1/final.png", width: 60%),
  caption: [
    Render of the geometry of the water molecule,
    optimized using MACE-ICE13-1.
  ],
)

=== Geometry optimization

To rapidly put to the test the calculators, geometrical values of the relaxed configuration have been computed and compared with references in the literature.
The bond lengths and the bond angle of the water molecule are known with remarkable accuracy from the vibration-rotation spectra of normal and isotopic water vapour. @eisenbergStructurePropertiesWater2005[§1.1 (c)]

The first value concerns the characteristic $#ce("HOH")$ *bend angle* of the molecule;
the accurate description of this physical value is a required test to ensure the validity of the models.
As can be seen in @table:hoh-angle, the MACE-MP-0 large model most accurately describes the bend angle of the molecule, while the small model is the most distant from reference @eisenbergStructurePropertiesWater2005[Table 1.3, data from Benedict _et al._ (1956) and Herzberg (1950)] @PhysicalChemistryWater2020.

The second value is the $#ce("OH")$ *bond length*.
@table:oh-bond-length shows that MACE-ICE13-1 gives the most accurate description according to reference @eisenbergStructurePropertiesWater2005[Table 1.2, data from Benedict _et al._ (1956)] @PhysicalChemistryWater2020, and MACE-MP-0 small again is less accurate than the other models.

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
    width: 70%,
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
The optimiser accepts two important input parameters. The first if _fmax_, the force threshold, defined in $"eV" dot angstrom^(-1)$.
The convergence criterion is that the force on all individual atoms should be less than fmax:
#footnote[https://wiki.fysik.dtu.dk/ase/ase/optimize.html]
$
  max_a |arrow(F)_a| < f_"max"
$

For the present purposes, a value of $f_"max"$ between $10^(-4)$ and $10^(-8) "eV/"angstrom$ yielded satisfactory results for the different calculator models tested.

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

=== Molecular vibrations and assessment of the dynamical stability <sec:molecule-vibrations>

The nuclei of molecules, far from occupying fixed positions with respect to each other, are in continual state of vibration.
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
$ <eq:molecule-energy-above-vibrationless-equilibrium-state>
where the sums are over normal modes.
The $omega$s in this equation are often called the _harmonic frequencies_;
they are the frequencies with which the molecule would vibrate if its vibrations were perfectly harmonic.
The $x$s are the _anharmonic constants_ and describe the effect on the vibrational frequencies of the departure from purely harmonic form of the vibrations.
@table:molecule-omega and @table:molecule-omega-errors
contain the harmonic frequencies for $#ce("H2O")$ and the errors for each model compared with reference data determined by Benedict _et al._ (1956) @eisenbergStructurePropertiesWater2005[Table 1.4].
The value of the anharmonic frequencies is reported for completeness in @table:molecule-vibrational-constants-anharmonic.

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

#figure(
  table(
    columns: 2,
    table.header(
      [Molecule],
      $ce("H2O")$,
    ),

    $x_11$, $-42.576$,
    $x_22$, $-16.813$,
    $x_33$, $-47.566$,
    $x_12$, $-15.933$,
    $x_13$, $-165.824$,
    $x_23$, $-20.332$,
  ),
  caption: [
    Anharmonic vibrational constants of $ce("H2O")$ for @eq:molecule-energy-above-vibrationless-equilibrium-state, taken from @eisenbergStructurePropertiesWater2005[Table 1.4].
    Units are in $"cm"^(-1)$.
  ],
) <table:molecule-vibrational-constants-anharmonic>

Studying the vibrational properties of the geometry obtained at the end of the optimization procedure also allows us to assess if the final geometry is a stable or unstable configuration.
The vibrational modes are calculated from a finite difference approximation of the Hessian matrix, displacing atoms according to a parameter named `delta`, measured in $angstrom$.
#footnote[https://wiki.fysik.dtu.dk/ase/ase/vibrations/modes.html]

The following figures represent the frequencies of the vibration modes of the optimized water molecule obtained with @ase plotted against the displacement parameter `delta`, using different calculator models and values of the `fmax` parameter.
Frequencies are indexed in ascending order, and *imaginary frequencies*, representing unstable configurations, are shown as negative values in the graphs.
Inclusion of dispersion contributions in calculators leads to minimal differences in the
converged configurations.

// logseq://graph/softseq?block-id=6651c91e-359a-4b31-9a4f-9dbb3ff5da83
The stability of configurations is dependent in a discriminant way on the `delta` and `fmax` parameters.
In particular, we seek simulation parameters that bring the energies and the spatial frequencies of the vibrations other than the first three, as close to zero as possible.
The analysis confirms that, to achieve the result with the best numerical accuracy possible, the value of $f_"max" = 10^(-8) "eV/"angstrom$ is the best.
However, such a value of the force threshold could be considered extreme, especially compared to typical values of $f_"max" = 0.01 "eV/"angstrom$ commonly used in @dft calculations.

With the same spirit, the displacement `delta` shall be smaller than $10^(-4) angstrom$ to ensure convergence of calculations of frequencies within precision of $10^(-2) " cm"^(-1)$,
and smaller than $10^(-3) angstrom$ for a convergence of the frequencies within $10^(-1) " cm"^(-1)$, which is well within the acceptable range.
Again, this is a extremely small value.
The default $delta = 0.01 angstrom$ in @ase is a pretty conservative displacement by solid-state standards, and the tight value obtained here could be explained by the fact that the $ce("H -")$ stretch bonds encountered here are very strong and smaller values of delta are better for the system considered in our calculations.

When seeking the force threshold, we should ultimately be guided by the following considerations:

- the force tolerance of the optimizer may need to be very tight;
- even then, the finite difference scheme is unlikely to get the trans/rotational modes to be _exactly_ zero;
  usually one aims for "small enough" values, which may be negative;
- a "noisy" potential energy surface is more likely to have problems with consistency between a local minimum and finite-difference vibrations;
  with #glspl("mlp") this can be an issue, especially if using single-precision GPU evaluation;
- on the other hand, a well implemented @mlp can be regularized to a smooth surface and may behave "better" than the underlying @dft calculations under tighter thresholds.


The MACE models produce smooth energy profiles @batatiaMACEHigherOrder2022[§5.3.2], which lends support to the hypothesis that smaller thresholds yield better results.

Relaxing the requirement for the obtainment of numerically exact zero values for the energies and frequencies of the vibrations other than the three normal modes, allows for the selection of less stringent values for $f_"max"$ and $delta$, depending on the desired accuracy.
A summary of the results is detailed in @table:vibrations-range-of-values.
In conclusion, a value of $f_"max" = 10^(-2) "eV/"angstrom$ guarantees a convergence of the energy of vibrations within $6 "meV"$, and a value of $f_"max" = 10^(-3) "eV/"angstrom$ a convergence within $1.8 "meV"$, and a value of $f_"max" = 10^(-4) "eV/"angstrom$ a convergence within $0.6 "meV"$.

Below is also available a grid of plots comparing the convergence of frequencies using different models, force thresholds and displacements.
Adding the dispersion correction to the models in this step produces negligible differences in the results, so their graphs are omitted for brevity.

#figure(
  table(
    columns: 5,
    align: horizon,
    table.header($f_"max"$, $10^(-1)$, $10^(-2)$, $10^(-3)$, $10^(-4)$),
    [small],
    [$17.7i$ \ $delta = 10^(-3)$],
    [$2 ÷ 6$ \ $delta = 10^(-2)$],
    [$0.4i$ \ $delta = 10^(-3)$],
    [$0.4i$ \ $delta = 10^(-3)$],

    [medium],
    [$1.4i ÷ 3.8$ \ $delta = 10^(-3)$],
    [$3.8 ÷ 5i$ \ $delta = 10^(-2)$],
    [$0.3 ÷ 1.1$ \ $delta = 10^(-3)$],
    [$0.2 ÷ 0.4$ \ $delta = 10^(-3)$],

    [large],
    [$3.4 ÷ 11.1$ \ $delta = 10^(-2)$],
    [$2÷4$ \ $delta = 10^(-2)$],
    [$0.2 ÷ 1.8$ \ $delta = 10^(-2)$],
    [$0.4i ÷ 0.6$ \ $delta = 10^(-2)$],
  ),
  caption: [
    Table showing the range of values obtained for the energies of vibrations outside the first three normal modes, expressed in $"meV"$, for each MACE-MP-0 model on the rows, and for each chosen $f_"max"$ on the columns, expressed in $"eV/"angstrom$.
    Inside each cell, on the second row is also shown the corresponding least restrictive value of the displacement $delta$, expressed in $angstrom$, that guarantees convergence of results for each configuration.
  ],
) <table:vibrations-range-of-values>

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

==== Zero-point vibrational energy <sec:molecule-zpe>

@eq:molecule-energy-above-vibrationless-equilibrium-state also yields an expression for the zero-point energy of vibration @eisenbergStructurePropertiesWater2005[Eq. 1.4]:
$
  "ZPE" &:= G(0,0,0) \
  &= 1 / 2 (omega_1 + omega_2 + omega_3) +
  1 / 4 (x_11 + x_22 + x_33 + x_12 + x_13 + x_23).
$ <eq:zpe>
When the reference harmonic and anharmonic constants of @table:molecule-omega and @table:molecule-vibrational-constants-anharmonic are inserted in this equation, the zero-point energy of $ce("H2O")$ is found to be $4634.32 " cm"^(-1)$, or $0.575 "eV"$.
It has already been pointed out that the @zpe cannot be expressed in terms of the $omega_i$ only, and this implies that different scaling factors must be used for obtaining improved #glspl("zpe") or fundamental vibrations from harmonic frequencies, by comparing to experimental values. @baroneVibrationalZeropointEnergies2004[§II] @grevConcerningZeropointVibrational1991
Conscious of the fact that our calculations are done in the harmonic approximation#footnote[#gls("ase", long: false) Vibrations class computes ZPE using @eq:zpe-harmonic; see https://wiki.fysik.dtu.dk/ase/_modules/ase/vibrations/data.html#VibrationsData.get_zero_point_energy], we ought not compare our results with the full @zpe, but with the following harmonic component: @baroneVibrationalZeropointEnergies2004[Eq. 5]
$
  "ZPE"_H := 1 / 2 sum_i omega_i.
$ <eq:zpe-harmonic>
The experimental value for this quantity is $"ZPE"_H^"Expt." = 56.4 "kJ/mol" approx 0.585 "eV"$. @baroneVibrationalZeropointEnergies2004[Table 1]
The harmonic @zpe computed using #glspl("mlp") is compared with this value in @table:zpe.
The MACE-ICE13-1 model shows the best agreement, with a discrepancy of $0.02 "eV"$, with the MACE-MP-0 small model following righ after.

#let zero_point_energies_table = csv("simulazioni/02_water/01_molecule/zero_point_energies.csv")
#figure(
  table(
    columns: zero_point_energies_table.first().len(),
    // table.header(..zero_point_energies_table.first()),
    table.header([Model], [ $"ZPE"_H ("eV")$], [Discrepancy (eV)]),
    ..zero_point_energies_table.slice(1).flatten(),
  ),
  caption: [
    Zero-point energies in the harmonic approximation for the @mlp calculators and reference experimental value from @baroneVibrationalZeropointEnergies2004[Table 1].
    The last column describes the difference between calculated and experimental value.
    small, medium and large models refer to MACE-MP-0.
  ],
) <table:zpe>

The overall results so far indicate that the medium, large MACE-MP-0, and MACE-ICE13-1 approximate better the properties of the water molecule,
while the small MACE-MP-0 model is the worst of them in this regard, except for the case of the @zpe calculation.

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

Analogously to the procedure adopted for the analysis of vibrations of the water monomer in @sec:molecule-vibrations,
normal modes of the water dimer were analyzed to assess the dynamical stability of the optimized geometry.
A displacement of $10^(-2) angstrom$ or less is appropriate for the MACE-ICE13-1 and MACE-MP-0 large models to converge to stable configurations within a threshold tolerance in frequencies of about $1"cm"^(-1)$, while for MACE-MP-0 small and medium a displacement of $10^(-3) angstrom$ or less is required to respect the same threshold tolerance in frequencies.

#large_figure(
  grid(
    columns: 2,
    image("simulazioni/02_water/02_dimer/01_optimize/Grafici/small.png"),
    image("simulazioni/02_water/02_dimer/01_optimize/Grafici/medium.png"),

    image("simulazioni/02_water/02_dimer/01_optimize/Grafici/large.png"),
    image("simulazioni/02_water/02_dimer/01_optimize/Grafici/MACE-ICE13-1.png"),
  ),
  caption: [
    Values of the vibration frequencies of the water dimer after structure optimization, plotted against the displacement used for the vibrations analysis.
    Negative values on the y-axis represent imaginary frequencies.
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
MACE-ICE13-1 demonstrates the best overall adherence to the prediction of harmonic frequencies.

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

=== Zero Point Energy

As it was done in @sec:molecule-zpe, we make use of @eq:zpe-harmonic to get the @zpe composed of the harmonic frequencies.
@table:dimer-zpe shows the results obtained with our @mlp models and compares them with a reference value from the literature.
The ranking of accuracy of the models sees as before MACE-ICE13-1 as closest to the ground truth, with MACE-MP-0 small following right after.

#let dimer_zpe = csv("simulazioni/02_water/02_dimer/01_optimize/Analisi/zpe.csv")
#figure(
  table(
    columns: dimer_zpe.first().len(),
    table.header(..dimer_zpe.first()),
    ..dimer_zpe.slice(1).flatten()
  ),
  caption: [
    ZPE in units of eV of the water dimer in the harmonic approximation and comparison with reference @kalesckyLocalVibrationalModes2012 (named Kalescky in the table).
  ],
) <table:dimer-zpe>

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

= Results II: crystal structures <sec:results-2>
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

Normal modes of vibration are calculated using the so-called *small displacement method*.
This method is increasingly more accurate with bigger and bigger supercells.

For details on the tools used for phonons calculations, see @sec:tools-phonons.

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
Frequencies are significantly higher than reference; the source of this discrepancy is not clear;
varying volumes of the cell were tested and did not change the result;
another source of the issue could be a built-in deviation of the calculator.
This issue has yet to be fully investigated at the time of writing.

Calculation of phonons dispersion along the band path was also performed using the PHON code @alfePHONProgramCalculate2009 for consistency of calculations.
The obtained results are in accordance with reference and previous calculations, maintaining the charachteristic over-estimation of the phonon frequencies we saw before.

#figure(
  image(
    "simulazioni/02_water/04_crystal_phonons/phon/12.PHON_MACE-ICE13-1_S3/bandplot.svg",
    width: 70%,
  ),
  caption: [
    Phonons bandstructure of ice Ih, computed with PHON using MACE-ICE13-1.
  ],
)

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

#figure(
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
  caption: [Execution times with phonopy. @togoFirstprinciplesPhononCalculations2023],
)

#figure(
  tablem(
    ignore-second-row: false,
    [
      |supercell|time|device|
      |2|1m 30s|cuda|
      |4|fail (out of memory)|cuda|
      |4|7h 32m| cpu|
    ],
  ),
  caption: [
    Execution times with Phonons by ASE.
    This setup achieves the slowest of timings, as ASE does not take into account symmetries.
  ],
)

#figure(
  table(
    columns: 4,
    [supercell], [forces time], [dispersions time], [device],
    [3], [3m 22s], [22s], [cuda],
  ),
  caption: [
    Execution times with PHON. @alfePHONProgramCalculate2009
  ],
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

== Molecular Dynamics
Constant NVT @md simulations with Langevin thermostat #footnote[https://wiki.fysik.dtu.dk/ase/ase/md.html#module-ase.md.langevin
] were performed under varying external conditions.
The thermostat couples the system to an external heat bath at a fixed temperature.

=== Radial Distribution Function

A fundamental property is the #gls("rdf", long: true), which describes how particle density varies as a function of distance from a reference particle.
@rdf is crucial for understanding structural properties in both solids and liquids.
This is a key property that needs to be reproduced accurately if the model is to be trusted for reliable predictions.
Any model that fails to capture the correct @rdf may lead to unreliable results for both structural and dynamic properties.

In this work, we utilize constant NVT simulations, with the system being thermally regulated using the Langevin thermostat.
The detailed simulation settings include the specification of the temperature of the heat bath, according to the reference data taken from experiment, the definition of the time step of the simulation and of the number of steps to simulate, related to the total physical time elapsed in the simulated system, which ensures thermal equilibrium throughout the simulation.
We demonstrate that the MACE-ICE13-1 model correctly reproduces the @rdf for liquid water within error, as shown in @fig:rdf and @fig:rdf-macemp0-maceice131.

The @rdf of the thermalized states is shown in @fig:rdf, compared with reference data @skinnerBenchmarkOxygenoxygenPairdistribution2013 from X-ray diffraction experiment.
To guarantee thermal equilibrium and total recombination of atom positions and bonds in the liquid, the simulated phisical time shall not be less than 100ps.
Constant NPT simulations should be more appropriate for the computation of physical properties, but they are missing at the present time.
Further analysis can also be made on the study of the diffusion coefficient and the density of the system.

#figure(
  image("simulazioni/02_water/05_md/Grafici/rdf_oo_mace-ice13-1_100ps_nbins=40.svg"),
  caption: [Radial distribution function of oxygens in liquid water.],
) <fig:rdf>

While the MACE-ICE13-1 model performs well in reproducing structural properties, further improvements could be made by enhancing the model's accuracy in predicting dynamic properties such as diffusion coefficients.

#large_figure(
  grid(
    columns: 2,
    image("simulazioni/02_water/05_md/Grafici/rdf_oo_mace-mp-0_NVT_T=297.15_t=5ps.svg"),
    image("simulazioni/02_water/05_md/Grafici/rdf_oo_mace-ice13-1_nbins=40.svg"),
  ),
  caption: [Comparison of the RDFs obtained from MD simulations of liquid water using MACE-MP-0 and MACE-ICE13-1.],
) <fig:rdf-macemp0-maceice131>

#figure(
  image("simulazioni/02_water/05_md/Grafici/temperature_NVT.png", width: 80%),
  caption: [
    Example of temperature trend through a MD simulation of water with MACE-ICE13-1.
  ],
)

= Tools <sec:tools>
== Hardware: Ibisco
Simulations were performed on the @ibisco cluster provided by the Federico II University. @WikiArchit_ib_enIbisco
// #box(
//   stroke: 2pt + red,
//   inset: 1mm,
//   [
The architecture of the hybrid cluster of the @ibisco Data Center can be represented as a set of multiple layers.
The lowest layer of the architecture consists of the hardware, characterized by calculation and storage nodes;
in the upper level the application level, which allows users to submit their tasks.
The intermediate level of the architecture consists of the set of CUDA and MPI libraries which are capable of making the two levels communicate with each other.

=== The hardware level
The cluster comprises 36 nodes and 2 switches, placed in 4 racks of the Data Center.
They perform two functions: calculation and storage.
To support the calculation there are 128 GPUs, distributed among 32 nodes (4 GPUs per node).
To support storage, 320 TB are available distributed among 4 nodes (80 TB per node).
To ensure access to resources and low-latency broadband communication between nodes, the InfiniBand technology is used to provide a high-performance network.

=== The Compute Node Architecture
The cluster compute nodes are 32 Dell C4140s, each equipped with 4 NVIDIA V100 GPUs, 2 Ethernet ports at 10Gb/s each, 2 InfiniBand ports at 100Gb/s each, 2 Intel Gen 2 Xeon Gold CPUs, and 2 SATA 480 GB SSDs.
Each node is also equipped with 22 64 GB RAM memory modules, overall 1.375 TiB.
Each GPU is equipped with 32 GB RAM memory.
The nodes are divided into 3 differently sized sub-clusters.
//   ],
// )

== Framework: Atomic Simulation Environment

The #gls("ase", long: true) @larsenAtomicSimulationEnvironment2017 is a set of tools and Python modules for setting up, manipulating, running, visualizing and analyzing atomistic simulations.
@ase provides interfaces to different codes through `Calculators` which are used together with the central `Atoms` object and the many available algorithms in @ase.
The `Atoms` object contains the positions of the atoms and the properties of the cell.
Among the many applications of @ase, it allows for geometry optimization, static calculations of potential energy and forces, as well as molecular dynamics simulations.
Calculations with the MACE #glspl("mlp") were performed through the calculators interface of @ase.

=== Structure optimization

The optimization algorithms can be roughly divided into local optimization algorithms which find a nearby local minimum and global optimization algorithms that try to find the global minimum (a much harder task).
Most optimization algorithms available in @ase accept the same base parameters:
- The atoms object to relax
- A log file
- An attached trajectory object
Basic optimizers optimize only internal atomic positions.
Cell volume and shape can also be optimized in combination with Filter tools.

The local optimization algorithms available in @ase follow the convergence criterion, which asks that the force on all individual atoms should be less than $f_"max"$:
$
  max_a |va(F_a)| < f_"max".
$

The BFGS algorithm uses two quantities to decide where to move the atoms on each step:
- the forces on each atom, as returned by the associated calculator;
- the Hessian matrix, i.e. the matrix of second derivatives $pdv(E, x_i, x_j)$ of the total energy with respect to nuclear coordinates.
If the atoms are close to the minimum, such that the potential energy surface is locally quadratic, the Hessian and forces accurately determine the required step to reach the optimal structure.
The Hessian is very expensive to calculate a priori, so instead the algorithm estimates it by means of an initial guess which is adjusted along the way depending on the information obtained on each step of the structure optimization.

== Calculator: MACE <sec:mace>

MACE
@Batatia2022mace
@Batatia2022Design
is an equivariant message-passing graph tensor network where each layer encodes many-body information of atomic geometry.
At each layer, many-body messages are formed using a linear combination of a tensor product basis. @batatiaDesignSpaceEquivariant2022 @darbyTensorReducedAtomicDensity2023
This is constructed by taking tensor products of a sum of two-body permutation-invariant polynomials, expanded in a spherical basis.
*The final output is the energy contribution of each atom to the total potential energy.*
The MACE architecture is implemented in PyTorch and employs the e3nn library.

In the following sections the technical details of the internals of MACE will be exposed.
A brief introduction to the topics that follow, #glspl("mpnn") and #glspl("gnn"), is detailed in @sec:gnn.

=== MPNN Interatomic Potentials <sec:mpnn>
#glspl("mpnn", long: true) are a type of #gls("gnn", long: true) that parametrizes a mapping from a labeled graph to a target space, either a graph or a vector space. @batatiaMACEHigherOrder2022
When applied to parametrize properties of atomistic structures (materials or molecules), the graph is embedded in 3-dimensional (3D) Euclidean space, where each node represents an atom, and edges connect nodes if the corresponding atoms are within a given distance of each other.
The state of each node $i$ in layer $t$ of the @mpnn is represented by a tuple
$
  sigma_i^((t)) = (arrow(r)_i, z_i, arrow(h)_i^((t))),
$
where $arrow(r)_i in RR^3$ is the position of atom $i$; $z_i$ is the chemical element; $arrow(h)_i^((t))$ are its learnable features.
A forward pass of the network consists of multiple _message construction, update_ and _readout_ steps.
During message construction, a message $arrow(m)_i^((t))$ is created for each node by pooling over its neighbours:
$
  arrow(m)_i^((t)) = plus.circle.big_(j in cal(N) (i)) M_t (
    sigma_i^((t)), sigma_j^((t))
  ),
$
where $M_t$ is a learnable message function and $plus.circle.big_(j in cal(N) (i))$ is a learnable, permutation invariant pooling operation over the neighbours of atom $i$ (e.g., a sum).
In the update step, the message $arrow(m)_i^((t))$ is transformed into new features
$
  arrow(h)_i^((t+1)) = U_t (sigma_i^((t)), arrow(m)_i^((t))),
$
where $U_t$ is a learnable update function.
After $T$ message construction and update steps, the learnable readouts functions $cal(R)_t$ map the node states $sigma_i^((t))$ to the target, in this case the site energy of atom $i$,
$
  E_i = sum_(t=1)^T cal(R)_t (sigma_i^((t))).
$

=== Equivariant Graph Neural Networks
In _equivariant_ #glspl("gnn"), internal features $arrow(h)_i^((t))$ transform in a specified way under some group action. @batatiaMACEHigherOrder2022
When modelling the potential energy of an atomic structure, the group of interest is $O(3)$, specifying rotations and reflections of the particles; translation invariance is trivially incorporated through the use of relative distances.
A @gnn is called $O(3)$ equivariant if it has internal features that transform under the rotation $Q in O(3)$ as
$
  arrow(h)_i^((t)) (Q dot (arrow(r)_1, dots, arrow(r)_N)) = D(Q) arrow(h)_i^((
    t
  )) (arrow(r)_1, dots, arrow(r)_N),
$
where $D^L(Q) in RR^((2L + 1) times (2L + 1))$ is a Wigner D-matrix of order $L$.
A feature labelled with $L=0$ describes an invariant scalar.
Features labelled with $L>0$ describe equivariant features, formally corresponding to equivariant vectors, matrices or higher order tensors.
The features of _invariant_ models, such as SchNet and DimeNet, transform according to $D(Q) = bb(1)$, the identity matrix.
Models such as NequIP, equivariant transformer, PaiNN, or SEGNNs, in addition to invariant scalars, employ equivariant internal features that transform like vectors or tensors.

=== The MACE Architecture

The MACE model follows the general framework of #glspl("mpnn") outlined in @sec:mpnn.
The key innovation is a new message construction mechanism.
The messages $arrow(m)_i^((t))$ are expanded in a hierarchical body order expansion,
$
  va(m_i^((t)))
  &= sum_j arrow(u)_1 (sigma_i^((t)); sigma_j^((t))) \
  &+ sum_(j_1, j_2) arrow(u)_2 (
    sigma_i^((t)); sigma_(j_1)^((t)); sigma_(j_2)^((t))
  )
  + dots
  + sum_(j_1, dots, j_nu) arrow(u)_nu (
    sigma_i^((t)); sigma_(j_1)^((t)); dots; sigma_(j_nu)^((t))
  ),
$ <eq:batatiaMACEHigherOrder2022-7>
where the $arrow(u)$ functions are learnable, the sums run over the neighbours of $i$, and $nu$ is a hyper-parameter corresponding to the maximum correlation order, the body order minus 1, of the message function with respect to the states.

==== Message Construction
At each iteration, MACE first embeds the edges using a learnable radial basis $R_(k l_1 l_2 l_3)^((t))$, a set of spherical harmonics $Y_(l_1)^(m_1)$, and a learnable embedding of the previous node features $h_(j, tilde(k) l_2 m_2)^((t))$ using weights $W_(k tilde(k) l_2)^((t))$. @batatiaMACEHigherOrder2022
The $A_i^((t))$ features are obtained by pooling over the neighbours $cal(N)(i)$ to obtain permutation invariant 2-body features whilst, crucially, retaining full directional information, and thus, full information about the atomic environment:
$
  A_(i, k l_3 m_3)^((t))
  = sum_(l_1 m_1, l_2 m_2) C_(l_1 m_1, l_2 m_2)^(l_3 m_3) dot \
  dot sum_(j in cal(N)(i)) R_(k l_1 l_2 l_3)^((t)) (r_(j i)) Y_(l_1)^(m_1) (
    hat(r)_(j i)
  )
  sum_(tilde(k)) W_(k tilde(k) l_2)^((t)) h_(j, tilde(k) l_2 m_2)^((t)),
$ <eq:batatiaMACEHigherOrder2022-8>
where $C_(l_1 m_1, l_2 m_2)^(l_3 m_3)$ are the standard Clebsh-Gordan coefficients ensuring that $A_(i, k l_3 m_3)^((t))$ maintain the correct equivariance; $r_(j i)$ is the (scalar) interatomic distance, and $hat(r)_(j i)$ is the corresponding unit vector.
$R_(k l_1 l_2 l_3)^((t))$ is obtained by feeding a set of radial features that embed the radial distance $r_(j i)$ using Bessel functions multiplied by a smooth polynomial cutoff to a multi-layer perceptron.
In the first layer, the node features $h_j^((t))$ correspond to the (invariant) chemical element $z_j$.
Therefore, @eq:batatiaMACEHigherOrder2022-8 can be further simplified:
$
  A_(i, k l_1 m_1)^((1))
  = sum_(j in cal(N)(i)) R_(k l_1)^((1)) (r_(j i)) Y_(l_1)^(m_1) (
    hat(r)_(j i)
  ) W_(k z_j)^((1)).
$ <eq:batatiaMACEHigherOrder2022-9>
This simplified operation is much cheaper, making the computational cost of the first layer low.

The _key_ operation of MACE is the efficient construction of higher order features from the $A_i^((t))$-features.
This is achieved by first forming tensor products of the features, and then symmetrising:
$
  B_(i, eta_nu k L M)^((t))
  = sum_(arrow(l m)) cal(C)_(eta_nu, arrow(l m))^(L M)
  product_(xi=1)^nu sum_(tilde(k)) w_(k tilde(k) l_xi)^((
    t
  )) A_(i, tilde(k) l_xi m_xi)^((t)),
  quad arrow(l m) = (l_1 m_1, dots, l_nu m_nu)
$ <eq:batatiaMACEHigherOrder2022-10>
where the coupling coefficients $cal(C)_(eta_nu)^(L M)$ corresponding to the generalized Clebsh-Gordan coefficients ensuring that $B_(i, eta_nu k L M)^((t))$ are $L$-equivariant, the weights $w_(k tilde(k) l_xi)^((t))$ are mixing the channels $(k)$ of $A_i^((t))$, and $nu$ is a given correlation order.
$cal(C)_(eta_nu, arrow(l m))^(L M)$ is very sparse and can be pre-computed such that @eq:batatiaMACEHigherOrder2022-10 can be evaluated efficiently.
The additional index $eta_nu$ simply enumerates all possible couplings of $l_1, dots, l_nu$ features that yield the selected equivariance specified by the $L$ index.
The $B_i^((t))$-features are constructed up to some maximum $nu$.
This variable in @eq:batatiaMACEHigherOrder2022-10 is the order of the tensor product, and hence, can be identified as the order of the many-body expansion terms in @eq:batatiaMACEHigherOrder2022-7.
The computationally expensive multi-dimensional sums over all triplets, quadruplets, etc..., are thus circumvented and absorbed into @eq:batatiaMACEHigherOrder2022-9 and @eq:batatiaMACEHigherOrder2022-8.

The message $m_i^((t))$ can now be written as a linear expansion
$
  m_(i,k L M)^((t))
  = sum_nu sum_(eta_nu) W_(z_i k L, eta_nu)^((t)) B_(i, eta_nu k L M)^((t)),
$ <eq:batatiaMACEHigherOrder2022-11>
where $W_(z_i k L, eta_nu)^((t))$ is a learnable weight matrix that depends on the chemical element $z_i$ of the receiving atom and message symmetry $L$.
Thus, MACE implicitly constructs each term $u$ in @eq:batatiaMACEHigherOrder2022-7 by a linear combination of $B_(i, eta_nu k L M)^((t))$ features of the corresponding body order.

Under mild conditions on the two-body bases $A_i^((t))$, the higher order features $B_(i, eta_nu k L M)^((t))$ can be interpreted as a _complete basis_ of many-body interactions, which can be computed at a cost comparable to pairwise interactions.
Because of this, the expansion @eq:batatiaMACEHigherOrder2022-11 is _systematic_.
It can in principle be converged to represent any smooth $(nu + 1)$-body equivariant mapping in the limit of infinitely many features. @dussonAtomicClusterExpansion2022

==== Update
In MACE, the update is a linear function of the message and the residual connection: @batatiaMACEHigherOrder2022
$
  h_(i, k L M)^((t + 1))
  = U_t^(k L) (sigma_i^((t)), m_i^((t)))
  = sum_(tilde(k)) W_(k L, tilde(k))^((t)) m_(i, tilde(k) L M)
  + sum_(tilde(k)) W_(z_i k L, tilde(k))^((t)) h_(i, tilde(k) L M)^((t)).
$

==== Readout
In the readout phase, the invariant part of the node features is mapped to a hierarchical decomposition of site energies via readout functions: @batatiaMACEHigherOrder2022
$
  E_i = E_i^((0)) + E_i^((1)) + dots + E_i^((T)), quad "where" \
  E_i^((t)) = cal(R)_t (h_i^((t)))
  = cases(
    sum_(tilde(k)) W_("readout", tilde(k))^((t)) h_(i, tilde(k) 00)^((t)) & "if" t < T,
    "MLP"_"readout"^((t)) ({h_(i, k 00)^((t))}_k) & "if" t = T
  )
$
The readout only depends on the invariant features $h_(i, k 00)^((t))$ to ensure that the site energy contributions $E_i^((t))$ are invariant as well.
To maintain body ordering, MACE uses linear readout functions for all layers except the last, where it uses a one-layer multi-layer perceptron.

=== Performance of MACE

It has been @kovacsEvaluationMACEForce2023 that shows that MACE generally outperforms alternatives for a wide range of systems, including liquid water.
In those simulations, the many-body equivariant MACE model is an improvement with respect to the three-body atom-centered symmetry function-based feed-forward neural network model (BPNN), the three-body invariant message passing model REANN and the two-body equivariant message passing model NequIP.

MACE-MP-0 is a general-purpose @ml model, trained on a public database of 150k inorganic crystals, that is capable of running stable molecular dynamics on molecules and materials. @batatiaFoundationModelAtomistic2023
The model can be applied out of the box and as a starting or "foundation model" for almost any atomistic system of interest and is thus a step towards democratising the revolution of @ml force fields by lowering the barriers to entry. @batatiaFoundationModelAtomistic2023

=== Fine-tuning a custom model <sec:fine-tuning>

In this thesis, we studied and followed the procedure to create a custom made MACE model.
Current developments are undergoing in the field, experimenting on the various techniques to fine-tune custom MACE models. @kaurDataefficientFinetuningFoundational2024
Fine-tuning a MACE model implies starting from a so called "foundation model" of MACE, e.g. MACE-MP-0, and execute further training of the model on new data.
That data is comprised of desired geometries (preferably under specific thermodynamic conditions), and attached potential energies, forces and stresses of the system, computed with the @pes of choice, such as the one of a specific @dft functional.
The fine-tuning thus is the procedure that allows MACE to accurately reproduce the charachteristic behaviour of the chosen @pes through a new training procedure, resulting in a fine-tuned @mlp model.
For the present purposes, we chose a toy model to represent ice Ih at a pressure of 1.01325 bar and a temperature of 100.0 K.

Fine-tuning a MACE model is composed of three steps, detailed as follows:

1. *Sample the phase space* to obtain representatives of the system under the thermodynamic conditions we are interesed in.
  This was obtained employing NPT dynamics as provided by @ase to generate $tilde 100 "ps"$ of dynamics, from which uncorrelated structures can be extracted.
  Of this total, 50 representatives were randomly chosen from the images after thermalization.
  See @fig-finetune-sample-temperature for the thermalization pattern.
  // TODO menziona scelta tra termostato di Berensden e Parrinello, con citazioni.
2. *Compute the reference* values for the energies, forces and stresses.
  For this step we executed single point self-consistent calculations on each representative of the samples,
  using a @paw revPBE-D3 @dft functional through @vasp.
  The output files are then converted, joined and shuffled to compose the training and test sets for the machine learning procedure.
3. *Fine tune* the foundation model on the new dataset; in our case MACE-MP-0 was chosen.
  #footnote([
    Reference procedures are available at https://github.com/ACEsuit/mace?tab=readme-ov-file#finetuning-foundation-models
  ])
  The training and test sets are the input data for the training of the neural network that underlies MACE.
  The training spans 2000 epochs, at the end of which the compiled model is obtained, along with statistics.
  Plots showing neural network loss and error on energy versus epochs are available in @fig-finetune-epochs.

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
    Loss (left) and MAE of the energy (right) over epochs plots for the fine-tuning of a new model on ice Ih, with MACE-MP-0 small as foundation model.
  ],
) <fig-finetune-epochs>

== Phonons calculations <sec:tools-phonons>

For phonons dispersion calculations, three main tools were tested:
- @ase built-in Phonons class. The current implementation of phonons calculation in @ase #footnote[https://wiki.fysik.dtu.dk/ase/ase/phonons.html] is outclassed by Phonopy, partly because the former does not take account of symmetries. #footnote[See https://gitlab.com/ase/ase-workshop-discussion/-/issues/7#note_245747917 and https://gitlab.com/ase/ase/-/issues/1235]
- Phonopy @togoFirstprinciplesPhononCalculations2023 is not native to @ase, so it needed a conversion interface for the atoms positions and forces, that was found online.#footnote[https://gitlab.com/drFaustroll/lavello/-/blob/01296fa0b4530f07cfb97a02358bdfa289f8199f/phonons.py]
- PHON @alfePHONProgramCalculate2009 is not native to @ase, and it also needed a conversion interface for the atoms positions and forces; that was not found, so a conversion interface between PHON and ASE was written for our specific purposes.#footnote[https://github.com/visika/tesi-magistrale/blob/32d5e30b73990daacb712e37872c22b89586f570/simulazioni/02_water/04_crystal_phonons/phon/00.template/forces.py]

= Conclusions

In this thesis work we learned how to use and implement in our computations the MACE @mlp calculator.
This allowed fast and accurate prototyping of different molecules configurations, phonon dispersions and ice polymorphs lattice energies, with orders of magnitude shorter execution times with respect to DFT methods.
It is noteworthy that test simulations could be executed effortlessly on a high-end laptop CPU, for a tight feedback loop between script creation and debugging.
The results obtained using MACE show a satisfiying agreement with the base models onto which it is trained, with the exception of the phonon dispersions, which show a slight over-estimation of the energies, while the shape of the graphs remain accurate.
Overall #glspl("mlp") demonstrate a convenient accuracy/cost ratio;
the combination of this characteristic with the ability to fine-tune models to better reproduce the behaviour of particular structures of one's interest, makes this approach very versatile, and a precious addition in the toolbox of material scientists.

= Acknowledgments
I would like to express my gratitude to Flaviano Della Pia for his invaluable assistance in learning the new tools and for providing me with insightful advice and practical tutorials. I am also indebted to Andrea Zen for his meticulous guidance and prompt advice. Finally, I would like to thank Dario Alfè for imparting me with the fundamentals of knowledge in this field and for steering my thesis in the right direction.

This work has been funded by project code PIR01_00011 “IBISCo”, PON 2014-2020, for all three entities (INFN, UNINA and CNR).

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

