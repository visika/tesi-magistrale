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
I cristalli molecolari rappresentano un'area di studio significativa nella scienza dei materiali,
per la loro importanza in campi come la farmaceutica, l'elettronica, le energie rinnovabili.

In questo lavoro, abbiamo usato dei potenziali machine learning (MLP) per modellare struttura e dinamica dei cristalli molecolari, prendendo come caso studio l'acqua.

L'acqua oltre a essere di rilievo per moltissime branche della scienza, presenta molte anomalie, dovute al bilanciamento tra legami a idrogeno e forze di dispersione, che la rendono difficile da modellizzare con approcci computazionali.

== Uso dei MLP

Bisogna affrontare il compromesso tra costo computazionale e accuratezza.

La mia tesi si concentra sull'utilizzo di Machine Learning Potentials (MLP) per modellare i cristalli molecolari, usando l'acqua come esempio rappresentativo.
L'applicazione di MLP cattura interazioni intermolecolari complesse con l'accuratezza degli approcci ab initio ma a un costo computazionale ridotto.

In questo lavoro abbiamo testato differenti modelli MLP per la stima di alcune proprietà dell'acqua, tra le quali energie di reticolo, dispersione dei fononi, e dinamica a temperatura finita.

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
La teoria del funzionale densità (DFT) è una formulazione del problema elettronico a molti corpi della meccanica quantistica, in termini della densità di stato fondamentale al posto che rispetto alla funzione d'onda dello stato fondamentale.

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

== Polimorfi delle strutture cristalline

#grid(
  columns: 2,
  [
    La corretta stima delle energie di reticolo, assoluta e relativa, è di grande interesse per l'acqua e per i cristalli molecolari in generale.

    I metodi Diffusion Monte Carlo (DMC) restituiscono quantità più accurate per le energie.

    Tra i funzionali DFT revPBE restituisce le proprietà di struttura e dinamiche più accurate del ghiaccio Ih.
  ],
  image("../thesis/imgs/dellapia2022_f1.jpeg", height: 80%),
)


// Energie di reticolo dei diversi polimorfi.
== Energia di reticolo assoluta

#grid(
  columns: 2,
  image(
    "../simulazioni/02_water/03_ICE13_lattice_energies/absolute_lattice_energy.svg",
    height: 70%,
  ),
  [
    La quantità fisica solitamente considerata per stabilire la stabilità di un cristallo è la sua energia di reticolo assoluta:
    $
      E_"lattice" := E_"crystal" - E_"gas"
    $
    L'errore assoluto medio (MAE) tra MACE-ICE13-1 e DMC è $0.90 "kJ/mol"$.
    I potenziali di scambio e correlazione sono ritenuti buoni se hanno un $"MAE" lt.tilde 4 "kJ/mol"$.
  ],
)

== Energia di reticolo relativa

#grid(
  columns: 2,
  image(
    "../simulazioni/02_water/03_ICE13_lattice_energies/relative_lattice_energy.svg",
    height: 70%,
  ),
  [
    C'è interesse a catturare la stabilità relativa dei polimorfi del ghiaccio, cioè rispetto a una data fase cristallina invece che rispetto alla fase gassosa.
    Calcoliamo l'energia di reticolo relativa rispetto al ghiaccio esagonale Ih:
    $
      Delta E_"lattice"^x := E_"lattice"^x - E_"lattice"^"Ih"
    $
    che è indipendente dalla configurazione del monomero nella fase gassosa.
  ],
)

== Dispersione dei fononi

#let phonons_height = 230pt

#grid(
  columns: 2,
  gutter: 10pt,
  image(
    "../simulazioni/02_water/04_crystal_phonons/phon/16.MACE_geometrie_Flaviano/combined.svg",
    height: phonons_height,
  ),
  image(
    "../simulazioni/02_water/04_crystal_phonons/phon/15.revPBED3/combined.svg",
    height: phonons_height,
  ),
)

#grid(
  columns: (1fr, 1fr),
  align: center,
  row-gutter: 10pt,
  [MACE-ICE13-1], [revPBE-D3],
  [$approx 1÷4$ min], [$gt.tilde 5$h],
)

== Dinamica molecolare

Infine ci si è occupati di simulare l'acqua nello stato liquido.

#image("../simulazioni/02_water/05_md/Grafici/rdf_oo_mace-ice13-1_100ps_nbins=40.svg")

= Conclusioni
// Nelle conclusioni riassumere ciò che si è fatto;
// nelle presentazioni, spesso di chi è fuori dal nostro settore,
// non si capisce cosa si è fatto noi e cosa si è fatto gli altri.
// Riassumere in una frase ciò che si è imparato durante questa esperienza

== Riepilogo dei risultati

== Limiti e sviluppi futuri
- È necessario avere grandi quantità di dati in partenza per addestrare nuovi modelli
- Manca ancora il supporto alle simulazioni su GPU multiple

== Potenziali applicazioni

#focus-slide[
  Grazie per l'attenzione!
]

#show: appendix

= Appendice

== Atomic Simulation Environment (ASE) <touying:hidden>

L'Atomic Simulation Environment (ASE) è un insieme di strumenti e moduli Python per impostare, manipolare, eseguire, visualizzare e analizzare simulazioni atomistiche. Il codice è liberamente disponibile sotto licenza GNU LGPL.

#slide[
  ASE fornisce interfacce a diversi codici attraverso i Calculator, che vengono utilizzati insieme all'oggetto centrale Atoms e ai numerosi algoritmi disponibili in ASE.

  #image("imgs/ase_calculators.png")
]

#slide()[
  #set text(size: 18pt)
  ```python
  # Example: structure optimization of water molecule
  from ase import Atoms
  from ase.optimize import BFGS
  from mace.calculators import mace_mp
  from ase.io import write
  h2o = Atoms("HOH", positions=[(0, 0, 0), (0, 1, 0), (0, 0, 1)])
  h2o.calc = mace_mp(model="medium", dispersion=False, device="cuda")

  opt = BFGS(h2o, trajectory="optimization.traj")
  opt.run(fmax=0.02, steps=1000)
  # BFGS:   0  19:10:49    -31.435229     2.2691
  # BFGS:   1  19:10:50    -31.490773     0.3740
  # BFGS:   2  19:10:50    -31.492791     0.0630
  # BFGS:   3  19:10:51    -31.492848     0.0023
  write('H2O.xyz', h2o)
  h2o.get_potential_energy()
  # -31.492847800329216
  ```
]

== Costruzione e passaggio del messaggio nelle GNN <touying:hidden>

== Fine tuning di MACE