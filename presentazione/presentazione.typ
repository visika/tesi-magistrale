#import "@preview/touying:0.5.2": *
#import themes.university: *

// https://touying-typ.github.io/docs/themes/university

#show: university-theme.with(
  aspect-ratio: "16-9",
  config-info(
    title: [Properties of molecular crystals using machine learning potentials],
    author: [Mariano Mollo],
    date: [17 ottobre 2024],
    institution: [
      Università degli Studi di Napoli Federico II \
      Dipartimento di Fisica Ettore Pancini
    ],
  ),
)

#title-slide(
  title: text(36pt)[Properties of molecular crystals using machine learning potentials],
  author: (
    [Candidato:\ Mariano Mollo],
    [Relatori:\ Prof. Dario Alfè,\ Prof. Andrea Zen],
  ),
  logo: place(
    top + right,
    dx: 30pt,
    dy: -30pt,
    image("../thesis/imgs/University_Federico_II_Logo.svg", width: 80pt),
  ),
)

== Panoramica <touying:hidden>

#components.adaptive-columns(outline(title: none, indent: 1em))

= Introduzione

== Importanza dei cristalli molecolari

== Uso dei MLP

= Teoria

== DFT

== MLP

= Risultati

== Struttura dell'acqua e dimero

== Strutture cristalline e liquido

= Conclusioni

== Riepilogo dei risultati

== Limiti e sviluppi futuri

== Potenziali applicazioni

#focus-slide[
  Grazie per l'attenzione!
]

#show: appendix

= Appendice
