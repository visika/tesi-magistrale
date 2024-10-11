#import "@preview/touying:0.5.2": *
#import themes.university: *
#import "@preview/fletcher:0.4.5" as fletcher: diagram, node, edge
#import "@preview/pinit:0.2.0": *
#import "@preview/physica:0.9.3": *

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

#show figure.caption: it => text(fill: gray.darken(40%), size: 18pt, it.body)

#let pinit-highlight-equation-from(
  height: 2em,
  pos: bottom,
  fill: rgb(0, 180, 255),
  highlight-pins,
  point-pin,
  highlight-dy: -0.6em,
  body,
) = {
  pinit-highlight(
    ..highlight-pins,
    dy: highlight-dy,
    fill: rgb(..fill.components().slice(0, -1), 40),
  )
  pinit-point-from(
    fill: fill,
    pin-dx: 0em,
    pin-dy: if pos == bottom {
      0.8em
    } else {
      -0.6em
    },
    body-dx: 0pt,
    body-dy: if pos == bottom {
      -1.7em
    } else {
      -1.6em
    },
    offset-dx: 0em,
    offset-dy: if pos == bottom {
      0.8em + height
    } else {
      -0.6em - height
    },
    point-pin,
    rect(
      inset: 0.5em,
      stroke: (bottom: 0.12em + fill),
      {
        set text(fill: fill)
        body
      },
    ),
  )
}

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

// == Panoramica <touying:hidden>

// #components.adaptive-columns(outline(title: none, indent: 1em))

// = Introduzione

== Il ghiaccio e i suoi polimorfi

#grid(
  columns: 2,
  align: horizon,
  [
    - Classe dei cristalli molecolari
    - Importanza di studiare le diverse fasi cristalline del ghiaccio
    - Difficoltà: caratterizzare le componenti di interazione, Coulomb, van der Waals
    - Spesso i polimorfi differiscono solo per piccole differenze di energia ($approx 1 "kJ/mol"$)
    - Per caratterizzare il diagramma T-P, sono importanti le proprietà dinamiche.
  ],
  [
    #image("../thesis/imgs/dellapia2022_f1.jpeg", height: 80%)
  ],
)

#place(bottom + right)[
  Della Pia et al. 2022
]

== Potenziali tradizionali e machine learning

#slide[
  #grid(
    columns: (1fr, 1fr),
    [
      - Compromesso tra costo e accuratezza (DMC, CC, DFT, FF)
      - Esempio: trovare il cristallo più stabile tramite simulazioni di MD seguendo la PES, che ha a che fare con le posizioni degli ioni
      - Uso dei potenziali ML per modellare la superficie di energia potenziale (PES) racchiudendo nell'addestramento il ruolo degli elettroni
    ],
    [
      #uncover("1")[
        // #place(horizon)[
        #image("imgs/Interpolation_Data.svg")
        // ]
      ]
      #uncover("2")[
        #place(horizon, image("imgs/Interpolation_example_polynomial.svg"))
      ]
      #uncover("3")[
        #place(horizon, image("imgs/Potential_Energy_Surface_for_Water.png"))
        #place(
          bottom + right,
          dy: 2.5cm,
          text(17pt)[AimNature, CC BY-SA 3.0, via Wikimedia Commons],
        )
      ]
    ],
  )
]

== L'architettura di MACE: GNN

#slide[
  #image("../thesis/imgs/distill.pub.gnn-intro-caffeine-adiacency-graph.png")
]

== Modelli di MACE

#slide[
  Ci sono tre modelli base di MACE, che variano in base alla quantità di parametri che costituiscono il modello, e uno specializzato fornito dall'esterno per la mia tesi:

  - MACE-MP-0#footnote[Batatia et al. 2023, A foundation model for atomistic materials chemistry] (small, medium, large)
  - MACE-ICE13-1 (adattato con fine-tuning in particolare al ghiaccio, per gentile concessione di Flaviano Della Pia.)

  #place(right, image("imgs/flaviano.jpg", height: 40%))
]

= Risultati
Per l'implementazione pratica dei modelli MLP, si è selezionata l'acqua come cristallo molecolare da studiare.

L'acqua ha diverse proprietà peculiari, che mettono a dura prova i calcolatori che provino a simularne il comportamento.

Le anomalie che presenta sono dovute al bilanciamento tra legami a idrogeno e forze di dispersione, che la rendono difficile da modellizzare.

== Ottimizzazione della geometria
La procedura per ottimizzare la geometria degli atomi in un sistema richiede la corretta impostazione di alcuni parametri, come:

- il calcolatore e il tipo di dato numerico #uncover("2-")[(MACE-MP-0 medium, float64)]
#meanwhile
- la geometria iniziale e le proprietà della cella #uncover("3-")[($hat("HOH") = 90°$ nel vuoto)]
#meanwhile
- l'algoritmo di ottimizzazione #uncover("4-")[(BFGS)]
#meanwhile
- la soglia di convergenza delle forze residue su ciascun atomo $a$:
$
  max_a |arrow(F)_a| < f_"max" uncover("5-", = 10^(-4) "eV/"angstrom)
$
#meanwhile
- il numero massimo di passaggi dell'ottimizzatore #uncover("6-")[(1000)]

= Monomero #place(top + center,
  image("../simulazioni/02_water/01_molecule/MACE-ICE13-1/final.png")
)

== Ottimizzazione della geometria
#grid(
  columns: 2,
  image("../simulazioni/02_water/01_molecule/Grafici/angle_convergence_mace_mp_0_large.svg"),
  [
    Si sono ottimizzati la distanza del legame OH e l'angolo formato dalla molecola HOH, ottenendo valori coerenti con la letteratura, entro $0.01 angstrom$ e $0.5°$, rispettivamente.

    Si osserva ad esempio la soglia di convergenza nelle forze residue per il monomero,
    $f_"max" = 10^(-4) "eV/"angstrom$.

  ],
)

== Frequenze di vibrazione armoniche

#let crimson = rgb("#c00000")

#slide[
  #place(
    bottom + left,
    image("../thesis/imgs/eisenberg2005_fig1.1.gif", height: 50%),
  )
  $
    G(v_1, v_2, v_3)
    = sum_(i=1)^3 omega_i (v_i + 1 / 2)
    + #pin(1) sum_(i=1)^3 sum_(k >= i)^3 x_(i k) (v_i + 1 / 2) (v_k + 1 / 2) #pin(2)
  $

  $
    "ZPE" &:= G(0,0,0) \
    &= 1 / 2 (omega_1 + omega_2 + omega_3)
    + #pin(3) 1 / 4 (x_11 + x_22 + x_33 + x_12 + x_13 + x_23) #pin(4)
  $

  #pause

  #pinit-line(stroke: 3pt + crimson, start-dy: 30pt, end-dy: -30pt, 1, 2)

  #pause

  #pinit-line(stroke: 3pt + crimson, start-dy: 30pt, end-dy: -30pt, 3, 4)

  #pause

  #place(bottom + right)[
    I contributi anarmonici *non* sono trascurabili \ per lo studio delle proprietà dell'acqua cristallina.
  ]

]

#grid(
  columns: (2fr, 1fr),
  row-gutter: 10pt,
  [
    // #let molecule_omega_table = csv("../simulazioni/02_water/01_molecule/Analysis/omega.csv")
    // #table(
    //   columns: molecule_omega_table.first().len(),
    //   table.header(
    //   // ..molecule_omega_table.first(),
    //   [Model], $omega_1$, $omega_2$, $omega_3$
    // ),
    //   ..molecule_omega_table.slice(1).flatten()
    // )
    #table(
      columns: 5,
      table.header([Modello], $omega_1$, $omega_2$, $omega_3$, $"ZPE" ("eV")$),
      [Reference], [3832.17], [1648.47], [3942.53], [0.585],
      [small], [3722.10], [1485.10], [3841.80], [0.561],
      [medium], [3592.20], [1601.10], [3736.50], [0.554],
      [large], [3695.60], [1497.30], [3814.80], [0.558],
      [ICE13-1], [3702.30], [1607.80], [3807.50], [0.565],
    )
  ],
  image("../simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium fmax=1e-1.svg"),

  [
    Si sono poi calcolate le frequenze di vibrazione armoniche per verificare il loro accordo con i dati tabulati in referenza.
    Le frequenze sull'asse y negativo indicano instabilità della geometria.
    Il _Displacement_ è lo spostamento degli atomi usato dalla tecnica di stima delle frequenze.
  ],
  image("../simulazioni/02_water/01_molecule/Grafici/MACE-MP-0 medium fmax=1e-8.svg"),
)

= Dimero #place(top + center,
  image("../simulazioni/02_water/02_dimer/01_optimize/MACE-ICE13-1/final.png")
)

== Ottimizzazione della geometria

Si è ripetuta la stessa procedura di ottimizzazione della geometria e studio delle frequenze di vibrazione armoniche e ZPE per il dimero.

#image("../thesis/imgs/klopper-fig1.gif")

#let dimer_geometry_table = csv("../simulazioni/02_water/02_dimer/01_optimize/geometria.csv")
#place(
  bottom + right,
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
)

== Frequenze di vibrazione armoniche
#grid(
  columns: (1.5fr, 1fr),
  gutter: 10pt,
  figure(
    image("../simulazioni/02_water/02_dimer/01_optimize/Grafici/harmonic_frequencies_errors_barchart.svg"),
    caption: [Scarto delle frequenze rispetto alla referenza],
  ),
  [
    MACE-ICE13-1 risulta il modello che meglio stima in generale le frequenze armoniche del dimero d'acqua.

    #align(right)[
      #table(
        columns: 3,
        table.header(
          [Model],
          [ZPE (eV)],
          [Error],
        ),

        [Reference], [1.264], [0.000],
        [ICE13-1], [1.218], [-0.046],
        [small], [1.213], [-0.051],
        [medium], [1.210], [-0.054],
        [large], [1.200], [-0.064],
      )
    ]
  ],
)

== Binding energy
#grid(
  columns: 2,
  gutter: 10pt,
  image("../simulazioni/02_water/02_dimer/02_binding_energy/binding_energy.svg"),
  [
    La binding energy si calcola come

    $ Delta E_2 := E_2 - 2 E_1 $

    MACE-ICE13-1 rientra entro i valori tipici ottenuti con DFT:

    #align(right)[
      #image(
        "../thesis/imgs/mukhopadhyayWaterDimerII2018_fig5.png",
        height: 50%,
      )
    ]
  ],
)

= Strutture cristalline #place(top + center,
  image("../strutture/ICE13/Ih/Ih.png")
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
    C'è interesse a catturare la stabilità relativa dei polimorfi del ghiaccio, cioè rispetto a una data fase cristallina (ad es. ghiaccio Ih) invece che rispetto alla fase gassosa:
    $
      Delta E_"lattice"^x := E_"lattice"^x - E_"lattice"^"Ih"
    $
    che è indipendente dalla configurazione del monomero nella fase gassosa.
  ],
)

== Dispersione dei fononi


#slide[
  #image("../simulazioni/02_water/04_crystal_phonons/phonopy/Grafici/bandpath_gupta.svg")#pin("brillouin")
][
  Per il ghiaccio Ih è stata calcolata la dispersione dei modi normali del cristallo, che rappresentano oscillazioni collettive delle particelle nel sistema con frequenza $omega_(va(q),s)$, denominate _fononi_.

  #v(0.7cm)

  $
    #pin(1)D(arrow(q))#pin(2) :=
    1 / M sum_j e^(i#pin(5)arrow(q)#pin(6)dot arrow(r)_j^0)#pin(3)Phi(arrow(r)_j^0)#pin(4)
  $

  $
    omega^2#pin(7)arrow(epsilon)#pin(8)= D(arrow(q)) dot arrow(epsilon)
  $

  #pause

  #let fill = rgb(0, 180, 255)
  #pinit-highlight(
    1,
    2,
    dy: -0.65em,
    fill: rgb(..fill.components().slice(0, -1), 40),
  )

  #pinit-point-from(
    (1, 2),
    pin-dy: 15pt,
    body-dx: -150pt,
    body-dy: 10pt,
    offset-dx: 5pt,
    offset-dy: 40pt,
    fill: fill,
    text(fill: fill)[Matrice dinamica],
  )

  #pause

  #let fill = rgb(150, 90, 170)
  #pinit-highlight(
    3,
    4,
    dy: -0.7em,
    fill: rgb(..fill.components().slice(0, -1), 40),
  )

  #pinit-point-from(
    (3, 4),
    pin-dy: 20pt,
    offset-dx: 5pt,
    offset-dy: 110pt,
    body-dx: -200pt,
    fill: fill,
    text(
      fill: rgb(150, 90, 170),
      [Matrice delle costanti di forza],
    ),
  )

  #pause

  #pinit-highlight-equation-from(
    (5, 6),
    (5, 6),
    pos: top,
    height: 0.7em,
    fill: rgb(0, 255, 127),
    [Vettore d'onda],
  )

  #pause

  #let fill = rgb(255, 170, 0)
  #pinit-point-from(
    (7, 8),
    pin-dy: 10pt,
    pin-dx: -5pt,
    offset-dx: -100pt,
    body-dx: -190pt,
    fill: fill,
    text(fill: fill)[Vettore polarizzazione],
  )

  #pinit-highlight(
    7,
    8,
    dy: -0.65em,
    fill: rgb(..fill.components().slice(0, -1), 40),
  )

  #pause

  #pinit-point-from(
    "brillouin",
    pin-dy: -100pt,
    pin-dx: 100pt,
    offset-dy: -40pt,
    offset-dx: 50pt,
    body-dx: -50pt,
    [Per fare tutto ciò,\ si campionano i punti della\ zona di Brillouin.],
  )
]

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

#slide[
  #grid(
    columns: 2,
    image("../simulazioni/02_water/04_crystal_phonons/phonopy/heat_capacity_all_temps.svg"),
    [
      Lo studio delle frequenze armoniche ci permette di valutare alcune quantità fisiche di interesse, come ad esempio la capacità termica.

      Risulta evidente come il fine-tuning sulle caratteristiche del ghiaccio fornisca modelli più accurati nella caratterizzazione del materiale.
    ],
  )
]

= Liquido #place(top + center, dx: 10pt, dy: 20pt,
  image("../strutture/128_molecules/render.png")
)

== Dinamica molecolare
#slide[
  #grid(
    columns: (1fr, 1.2fr),
    [#image("../simulazioni/02_water/05_md/Grafici/temperature_NVT.png")
      Algoritmo di Verlet:],
    [
      La dinamica molecolare è un metodo che permette di campionare lo spazio delle fasi di un sistema isolato di $N$ particelle classiche interagenti, che obbediscono alle equazioni del moto di Newton:
      $
        va(f_i) = M dot.double(va(r))_i
        = - pdv(U({va(r)}), va(r_i))
      $
    ],
  )
  $
    arrow(r)_i (t + delta t)
    = arrow(r)_i (t)
    + arrow(v)_i (t) delta t
    + 1 / (2M) arrow(f)_i (t) (delta t)^2
    + cal(O)((delta t)^4)
  $
]

#slide[
  #image("../simulazioni/02_water/05_md/Grafici/rdf_oo_mace-ice13-1_100ps_nbins=40.svg")
][
  L'ultima task è lo studio del comportamento dell'acqua nello stato liquido.
  Si sono effettuate simulazioni di dinamica molecolare a NVT costanti con un termostato di Langevin, studiando la Radial Distribution Function (RDF), trovando un buon accordo con la referenza.
]

= Conclusioni
// Nelle conclusioni riassumere ciò che si è fatto;
// nelle presentazioni, spesso di chi è fuori dal nostro settore,
// non si capisce cosa si è fatto noi e cosa si è fatto gli altri.
// Riassumere in una frase ciò che si è imparato durante questa esperienza

== Riepilogo dei risultati

- Ho imparato come usare e implementare il calcolatore MLP MACE in simulazioni dei materiali, insieme agli strumenti a corredo (ASE, VASP, LAMMPS)
- Prototipazione veloce e accurata di diverse configurazioni molecolari, dispersione dei fononi, energie di reticolo di diversi polimorfi e dinamica molecolare, con tempi di esecuzione ordini di grandezza minori rispetto ai metodi DFT
- Potenziale con previsioni entro l'accuratezza chimica ($approx 4 "kJ/mol"$)
- Ho seguito la procedura di fine-tuning per creare un toy model specializzato di MACE
- Ho scritto script di integrazione tra i diversi strumenti, dove non esistevano (integrazione PHON-ASE)

== Limiti e sviluppi futuri
- È necessario avere grandi quantità di dati in partenza per addestrare nuovi modelli
- Per l'addestramento è necessario una quantità di risorse computazionali ingente (per l'utilizzo dei modelli no)
- Manca ancora il supporto alle simulazioni su GPU multiple
- Estensione degli studi ad altri cristalli molecolari e sistemi di interesse
- Democratizzazione dell'accesso alle simulazione accurata dei materiali

MACE risulta un'ottima aggiunta alla cassetta degli attrezzi nella scienza dei materiali.

== Potenziali applicazioni

#focus-slide[
  Grazie per l'attenzione!
]

#show: appendix

= Appendice

== Importanza dei cristalli molecolari

#slide[
  I cristalli molecolari rappresentano un'area di studio significativa nella scienza dei materiali,
  per la loro importanza in campi come la farmaceutica, l'elettronica, le energie rinnovabili.

  È importante il ruolo delle vibrazioni e delle interazioni van der Waals a molti corpi nelle proprietà coesive dei cristalli molecolari.
][
  #align(
    right,
    image("../thesis/imgs/karki2009_graphical_abstact.jpg", height: 50%),
  )
]

== Diversi potenziali per le simulazioni

#grid(
  columns: 2,
  gutter: 10pt,
  [
    #image("../thesis/imgs/images_large_cr0c00868_0001.jpeg")
  ],
  [
    - Compromesso tra costo computazionale e accuratezza.
    - Uso dei Machine Learning Potential (MLP) per modellare cristalli molecolari
    - Catturano interazioni intermolecolari complesse con l'accuratezza degli approcci ab initio ma a un costo computazionale ridotto.
  ],
)

In questo lavoro abbiamo testato differenti modelli MLP per la stima di alcune proprietà dell'acqua, tra le quali energie di reticolo, dispersione dei fononi, e dinamica a temperatura finita.

== Problematica e obiettivo

Ci sono difficoltà nei metodi tradizionali nel bilanciare l'accuratezza e il costo computazionale, come la Density Functional Theory (DFT), che pur essendo accurata, è costosa in termini di calcolo.

#grid(
  columns: 2,
  gutter: 20pt,
  image("../thesis/imgs/gilmerNeuralMessagePassing2017_Figure1.png"),
  [
    Vogliamo capire se è possibile migliorare la simulazione dei cristalli molecolari con accuratezza simile ai metodi ab initio ma a un costo computazionale molto più basso.
  ],
)

== DFT

La teoria del funzionale densità (DFT) è una formulazione in principio esatta del problema elettronico a molti corpi della meccanica quantistica, in termini della densità dello stato fondamentale degli elettroni, invece che della funziona d'onda di stato fondamentale.

Una delle sue numerose applicazioni è la determinazione di proprietà strutturali ed elettroniche dei materiali nella fisica dello stato solido.

L'uso di questi metodi richiede l'individuazione del corretto funzionale di scambio e correlazione.
È importante includere interazioni di van der Waals e di dispersione, allorquando non siano incluse nella teoria.

#place(bottom + right)[
  $
    hat(h)_"KS" [rho] (va(r)) psi_n (va(r)) =
    [-1 / 2 nabla^2 + v_"KS" [rho] (arrow(r))] psi_n (arrow(r))
    = epsilon_n psi_n (va(r))
  $
]

== Reti neurali a grafo

#slide[
  - Ciascun nodo rappresenta un atomo
  - Gli archi connettono i nodi se i corrispondenti atomi sono entro una data distanza tra loro
  - Lo stato di ciascun nodo $i$ nel layer $t$ è rappresentato da una tupla:

  #v(0.5cm)

  $
    sigma_i^((t)) = (
      #pin("r1")va(r_i)#pin("r2"),#pin("z1")z_i#pin("z2"),#pin("h1")va(h_i^((t)))#pin("h2")
    )
  $
][
  #image("../thesis/imgs/distill.pub.gnn-intro-caffeine-adiacency-graph.png")
  #image("../thesis/imgs/distill.pub.gnn-intro.fig1.png")

  #pause

  #pinit-highlight(
    "r1",
    "r2",
    dy: -0.65em,
  )

  #pinit-point-from(
    ("r1", "r2"),
    pin-dy: 15pt,
    pin-dx: 0pt,
    offset-dx: -15pt,
    body-dx: -200pt,
    [Posizione dell'atomo $i$],
  )

  #pause

  #pinit-highlight-equation-from(
    ("z1", "z2"),
    ("z1", "z2"),
    [Specie chimica],
  )

  #pause

  #let fill = rgb(150, 90, 170)
  #pinit-highlight(
    "h1",
    "h2",
    dy: -0.8em,
    dx: -0.1em,
    fill: rgb(..fill.components().slice(0, -1), 40),
  )

  #pinit-point-from(
    ("h1", "h2"),
    pin-dy: -20pt,
    pin-dx: 0pt,
    offset-dy: -60pt,
    offset-dx: -10pt,
    body-dx: -130pt,
    body-dy: -20pt,
    fill: fill,
    text(fill)[Caratteristiche \ apprendibili],
  )
]

#slide[
  #grid(
    columns: 2,
    row-gutter: 40pt,
    column-gutter: 20pt,
    [
      - Un passaggio tra layer consiste nella _costruzione_, _aggiornamento_ e _lettura_ di _messaggi_
    ],
    $
      va(m_i^((t))) = plus.circle.big_(j in cal(N) (i)) M_t (
        sigma_i^((t)), sigma_j^((t))
      )
    $,

    [
      - Nella fase di aggiornamento, il messaggio $va(m_i^((t)))$ è trasformato in nuove caratteristiche
    ],
    $
      arrow(h)_i^((t+1)) = U_t (sigma_i^((t)), arrow(m)_i^((t)))
    $,

    [
      - L'output finale è il contributo di energia di ciascun atomo all'energia potenziale
    ],
    $
      E_i = sum_(t=1)^T cal(R)_t (sigma_i^((t)))
    $,
  )

  Tutti gli operatori sono apprendibili ($M_t$, $plus.circle.big_(j in cal(N) (i))$, $U_t$, $cal(R)_t$)
]

== Proprietà di MACE

#slide[
  - È equivariante: le caratteristiche interne $va(h_i^((t)))$ si trasformano in modo preciso sotto l'azione di un gruppo, e.g. $O(3)$:
  $
    va(h_i^((t))) (Q dot (va(r_1), dots, va(r_N)))
    =pin("d1")D(Q)pin("d2")va(h_i^((t))) (va(r_1), dots, va(r_N))
  $

  - È message-passing, e costruisce messaggi secondo uno schema originale:

  $
    va(m_i^((t)))
    &= sum_j va(u_1) (sigma_i^((t)); sigma_j^((t))) \
    &+ sum_(j_1, j_2) va(u_2) (
      sigma_i^((t)); sigma_(j_1)^((t)); sigma_(j_2)^((t))
    )
    + dots
    + sum_(j_1, dots, j_nu)pin("u1")va(u_nu)pin("u2")(
      sigma_i^((t)); sigma_(j_1)^((t)); dots; sigma_(j_nu)^((t))
    )
  $

  #pinit-highlight-equation-from(
    ("u1", "u2"),
    ("u1", "u2"),
    pos: top,
    [Apprendibili],
  )

  #pinit-highlight-equation-from(
    ("d1", "d2"),
    ("d1", "d2"),
    pos: top,
    height: 20pt,
    highlight-dy: -0.7em,
    fill: rgb(150, 90, 170),
    [#h(3cm) Matrice D di Wigner],
  )
]

== Atomic Simulation Environment (ASE)

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

== Costruzione e passaggio del messaggio nelle GNN

== Fine tuning di MACE

#place(
  bottom + left,
  dx: -50pt,
  image(
    "../tutorial-fine-tuning/1.generate-training/run_2024-06-08/Figure_1_i_T.png",
    height: 55%,
  ),
)

#place(
  bottom + left,
  dx: 270pt,
  image("../tutorial-fine-tuning/analysis/loss_over_epochs.svg", height: 55%),
)

#place(
  bottom + left,
  dx: 470pt,
  image(
    "../tutorial-fine-tuning/analysis/mae_e_per_atom_over_epochs.svg",
    height: 55%,
  ),
)

#figure([
  #let my-node(..args) = node(
    corner-radius: 5pt,
    ..args,
  )

  #diagram(
    node-stroke: 1pt,
    my-node((0, 0), [Genera campioni \ da NPT MD]),
    edge(
    "-|>",
    [data_for_train.extxyz],
    bend: 45deg,
    // label-sep: 1em,
  ),
    my-node((1, 0), [Calcola riferimento \ da DFT]),
    edge(
      "-|>",
      align(
        center,
        [training_set.xyz \ test_set.xyz],
      ),
      bend: 45deg,
    ),
    my-node((2, 0), [Fine tune \ di MACE]),
    edge("-|>"),
    my-node((3, 0.5), [Modello]),
    edge((2, 0), (3, -0.5), "-|>"),
    my-node((3, -0.5), [Errori]),
  )
])

== Andamenti di tempi e memoria
#grid(
  columns: 2,
  image("../simulazioni/03_mace_memory_usage/execution_time_cpu_gpu.png"),
  image("../simulazioni/03_mace_memory_usage/execution_time_gpu.png"),

  image("../simulazioni/03_mace_memory_usage/max_memory_usage_cpu_gpu.png"),
  image("../simulazioni/03_mace_memory_usage/max_memory_usage_gpu.png"),
)
