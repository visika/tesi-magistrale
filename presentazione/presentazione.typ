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
Importanti nei campi della farmaceutica, dell'elettronica, delle energie rinnovabili.

== Uso dei MLP

La mia tesi si concentra sull'utilizzo di Machine Learning Potentials (MLP) per modellare i cristalli molecolari, usando l'acqua come esempio rappresentativo.

== Problematica e obiettivo

Ci sono difficoltà nei metodi tradizionali nel bilanciare l'accuratezza e il costo computazionale, come la Density Functional Theory (DFT), che pur essendo accurata, è costosa in termini di calcolo.

#grid(
  columns: 2,
  gutter: 20pt,
  image("../thesis/imgs/gilmerNeuralMessagePassing2017_Figure1.png"),
  [
    Vogliamo capire se è possibile migliorare la simulazione dei cristalli molecolari con accuratezza simile ai metodi ab initio ma a un costo computazionale molto più basso.

    // I MLP ci permettono di eseguire simulazioni di fisica con grande risparmio di tempo.
    // Ma a quale costo?

    // Investighiamo la questione confrontando i MLP con i modelli ab-initio sulla quale si basano.
  ],
)


= Teoria

== DFT
- revPBE-D3
- DMC

== MLP

#grid(
  columns: 2,
  gutter: 20pt,
  image("../thesis/imgs/distill.pub.gnn-intro-caffeine-adiacency-graph.png"),
)

I potenziali machine learning sfruttano le reti neurali per modellare le superfici di energia potenziale con alta precisione.

Per come sono strutturate, le reti neurali a grafo (GNN) sono molto versatili per questo scopo.

== MACE
È l'architettura specifica impiegata nel corso di questa tesi.

= Risultati
Applicazione di modelli MLP per predire proprietà dei cristalli molecolari.

Il nostro case study è l'acqua.
Per le sue proprietà peculiari, che mettono a dura prova i calcolatori che provino a simularne il comportamento.

== Struttura dell'acqua e dimero

== Frequenze di vibrazione armoniche

#image("../simulazioni/02_water/02_dimer/01_optimize/Grafici/harmonic_frequencies_errors_barchart.svg")

== Binding energy
#grid(
  columns: 2,
  image("../simulazioni/02_water/02_dimer/02_binding_energy/binding_energy.svg"),
  image("../thesis/imgs/mukhopadhyayWaterDimerII2018_fig5.png"),
)

== Strutture cristalline e liquido
Energie di reticolo dei diversi polimorfi.

#image("../simulazioni/02_water/03_ICE13_lattice_energies/absolute_lattice_energy.svg")

== Relative lattice energy

#image("../simulazioni/02_water/03_ICE13_lattice_energies/relative_lattice_energy.svg")

== Dispersione dei fononi

#grid(
  columns: 2,
  figure(
    image("../simulazioni/02_water/04_crystal_phonons/phon/16.MACE_geometrie_Flaviano/FREQ1.svg"),
    caption: [MACE-ICE13-1],
  ),
  figure(
    image("../simulazioni/02_water/04_crystal_phonons/phon/15.revPBED3/FREQ1.svg"),
    caption: [revPBE-D3],
  ),
)

== Dinamica molecolare

#image("../simulazioni/02_water/05_md/Grafici/rdf_oo_mace-ice13-1_100ps_nbins=40.svg")

= Conclusioni
// Nelle conclusioni riassumere ciò che si è fatto;
// nelle presentazioni, spesso di chi è fuori dal nostro settore,
// non si capisce cosa si è fatto noi e cosa si è fatto gli altri.
// Riassumere in una frase ciò che si è imparato durante questa esperienza

== Riepilogo dei risultati

== Limiti e sviluppi futuri
- Necessario sempre avere grandi quantità di dati in partenza per addestrare nuovi modelli

== Potenziali applicazioni

#focus-slide[
  Grazie per l'attenzione!
]

#show: appendix

= Appendice

== Costruzione e passaggio del messaggio nelle GNN <touying:hidden>
