#import "@preview/modern-unito-thesis:0.1.0": template

// Your acknowledgments (Ringraziamenti) go here
#let acknowledgments = [ 
  I would like to thank you for using my template and the team of typst for the great work they have done and the awesome tool they developed. Remember that it's best practice to thank the people you worked with before thanking your family and friends.
]

// Your abstract goes here
#let abstract = [
Molecular crystals play an important role in the field of materials science, particularly in drug development, electronics, and renewable energy sectors.

In this work we will study the properties of molecular crystals, using recently developed Machine Learning potentials to model their behaviour and characteristics. We will be primarily focusing on water as the initial subject, followed by a study of a selection of other molecular crystals.

Traditional approaches often grapple with the trade-off between computational expense and accuracy. The application of Machine Learning potentials captures complex intermolecular interactions with a significantly reduced computational cost compared to traditional ab-initio methods.

We will study the capabilities of trained Machine Learning potentials to accurately predict lattice energies, polymorphic behaviours, and response to external conditions like temperature and pressure. We will also study dynamic properties such as phonon spectra to complete the insight into the physical and chemical behaviours of molecular crystals.
]

#show: template.with(
  // Your title goes here
  title: "Properties of molecular crystals using machine learning potentials",

  // Change to the correct academic year, e.g. "2024/2025"
  academic-year: [2023/2024],

  // Change to the correct subtitle, i.e. "Tesi di Laurea Triennale",
  // "Master's Thesis", "PhD Thesis", etc.
  subtitle: "Master's Thesis",

  // Change to your name and matricola
  candidate: (
    name: "Mariano Mollo",
    matricola: "N94000618"
  ),

  // Change to your supervisor's name
  supervisor: (
    "Prof. Dario Alfè"
  ),

  // Add as many co-supervisors as you need or remove the entry
  // if none are needed
  co-supervisor: (
      "Prof. Andrea Zen",
  ),

  // Customize with your own school and degree
  affiliation: (
    university: "Università degli Studi di Napoli Federico II",
    school: "Scuola Politecnica e delle Scienze di Base",
    degree: "Corso di Laurea Magistrale in Fisica",
  ),

  // Change to "it" for the Italian template
  lang: "en",

  // University logo
  logo: image("imgs/University_Federico_II_Logo.svg", width: 40%),

  // Hayagriva bibliography is the default one, if you want to use a
  // BibTeX file, pass a .bib file instead (e.g. "works.bib")
  bibliography: bibliography("works.yml"),

  // See the `acknowledgments` and `abstract` variables above
  acknowledgments: acknowledgments,
  abstract: abstract,

  // Add as many keywords as you need, or remove the entry if none
  // are needed
  keywords: [Molecular Crystals, MACE, Molecular Dynamics]
)

// I suggest adding each chapter in a separate typst file under the
// `chapters` directory, and then importing them here.

#include "chapters/introduction.typ"

#include "chapters/example.typ"

#include "chapters/conclusions.typ"
