// PeakLocus workflow: inputs, five steps, one transcript, six region classes, example output.
#set page(width: 16cm, height: auto, margin: 0.8cm, fill: white)
#set text(font: "Libertinus Serif", size: 10pt, fill: rgb("#1d1d1b"))

#let ink = rgb("#1d1d1b")
#let muted = rgb("#5c5c58")
#let rule-c = rgb("#c9c9c2")
#let c-utr5 = rgb("#4f9a5a")
#let c-cds = rgb("#2f6f6a")
#let c-utr3 = rgb("#c27c2c")
#let c-exon = rgb("#4a7fb0")
#let c-intron = rgb("#9a9a94")
#let c-na = rgb("#d9d9d4")

#let heading(t) = text(size: 11pt, weight: "bold", fill: ink)[#t]

// Gene example: schematic coordinates 0 to 100 across the track.
#let tw = 14.4cm
#let X(v) = v * tw / 100
#let cx(i) = (i + 0.5) * tw / 5

#let seg(x0, x1, fill, label, dark: false) = place(
  dx: X(x0),
  dy: 0pt,
  rect(width: X(x1 - x0), height: 22pt, fill: fill, stroke: none, inset: 0pt)[
    #align(center + horizon)[#text(size: 8.5pt, weight: "bold", fill: if dark { white } else { ink })[#label]]
  ],
)

#let intron-line(x0, x1) = place(dx: X(x0), dy: 10pt, line(length: X(x1 - x0), stroke: 0.9pt + c-intron))

#let step(title, body) = pad(x: 5pt)[
  #box(width: 100%)[
  #align(center)[#text(size: 9pt, weight: "bold")[#title]]
  #v(1pt)
  #align(center)[#text(size: 8pt, fill: muted)[#body]]
  ]
]

#let chip(fill, name, body) = grid(
  columns: (14pt, 1fr),
  column-gutter: 6pt,
  rect(width: 14pt, height: 14pt, fill: fill, stroke: none),
  [#text(size: 9pt, weight: "bold")[#name] #linebreak() #text(size: 8.5pt, fill: muted)[#body]],
)

// Title
#text(size: 30pt, weight: "bold", fill: rgb("#1b3a5c"))[PeakLocus]
#v(-4pt)
#text(size: 12pt, fill: muted)[Transcript-aware genomic peak annotation]

#v(0.6cm)

// Inputs
#heading[Inputs]
#v(4pt)
#grid(
  columns: (1fr, 1fr, 1fr),
  column-gutter: 8pt,
  rect(width: 100%, stroke: 0.6pt + rule-c, inset: 7pt)[
    #text(size: 9pt, weight: "bold")[Peak intervals] #linebreak()
    #text(size: 8.5pt, fill: muted)[BED, narrowPeak, or tabular]
  ],
  rect(width: 100%, stroke: 0.6pt + rule-c, inset: 7pt)[
    #text(size: 9pt, weight: "bold")[Gene annotation] #linebreak()
    #text(size: 8.5pt, fill: muted)[GTF file]
  ],
  rect(width: 100%, stroke: 0.6pt + rule-c, inset: 7pt)[
    #text(size: 9pt, weight: "bold")[Strand, optional] #linebreak()
    #text(size: 8.5pt, fill: muted)[detected automatically if absent]
  ],
)

#v(0.6cm)

// Workflow
#heading[Workflow]
#v(6pt)
#box(width: tw, height: 18pt)[
  #for i in range(4) {
    place(dx: cx(i) + 11pt, dy: 8.5pt, line(length: cx(i + 1) - cx(i) - 24pt, stroke: 0.9pt + muted))
    place(dx: cx(i + 1) - 13pt, dy: 4pt, polygon(fill: muted, (0pt, 0pt), (0pt, 10pt), (7pt, 5pt)))
  }
  #for i in range(5) {
    place(dx: cx(i) - 9pt, dy: 0pt, circle(radius: 9pt, fill: c-cds, stroke: none)[
      #align(center + horizon)[#text(size: 9pt, weight: "bold", fill: white)[#(i + 1)]]
    ])
  }
]
#v(4pt)
#grid(
  columns: (1fr, 1fr, 1fr, 1fr, 1fr),
  step([Read peaks], [Detect columns and strand]),
  step([Load GTF], [Keep transcript features]),
  step([Normalize], [BED to GTF coordinates]),
  step([Map], [Overlap peaks with transcript regions]),
  step([Return], [One row per input peak]),
)

#v(0.6cm)

// Example transcript
#heading[Example transcript on the plus strand]
#v(8pt)
#box(width: tw, height: 36pt)[
  #place(dy: 2pt)[#seg(0, 12, c-utr5, [5′ UTR])]
  #place(dy: 2pt)[#seg(12, 50, c-cds, [CDS], dark: true)]
  #place(dy: 2pt)[#seg(50, 60, c-utr3, [3′ UTR])]
  #place(dy: 2pt)[#intron-line(60, 74)]
  #place(dy: 2pt)[#seg(74, 100, c-exon, [exon], dark: true)]
  #place(dx: X(62) - 16pt, dy: 30pt)[#text(size: 8pt, fill: muted)[intron]]
]

#v(0.6cm)

// Six classes
#heading[Region classes]
#v(6pt)
#grid(
  columns: (1fr, 1fr),
  column-gutter: 16pt,
  row-gutter: 8pt,
  chip(c-utr5, [5′ UTR], [5′ untranslated region, from the CDS bounds]),
  chip(c-utr3, [3′ UTR], [3′ untranslated region, from the CDS bounds]),
  chip(c-cds, [CDS], [Coding sequence]),
  chip(c-exon, [exon], [Exonic sequence not in a UTR or CDS, such as lncRNA]),
  chip(c-intron, [intron], [Between exons of one transcript]),
  chip(c-na, [N/A], [Overlaps no transcript feature]),
)

#v(0.6cm)

// Example output
#heading[Example output]
#v(6pt)
#set text(size: 8.5pt)
#table(
  columns: (auto, auto, auto, auto, auto, auto),
  stroke: (x, y) => if y == 0 { (bottom: 0.8pt + ink) } else { (bottom: 0.4pt + rule-c) },
  inset: 5pt,
  align: (left, left, left, center, left, left),
  table.header(
    [*Peak*], [*Gene*], [*Transcript*], [*Strand*], [*Location*], [*Gene type*],
  ),
  [chr1:104-110], [Gene1], [T1], [+], [5′UTR], [protein_coding],
  [chr1:129-140], [Gene1], [T1], [+], [CDS], [protein_coding],
  [chr1:249-260], [Gene1], [T1], [+], [intron], [protein_coding],
  [chr1:879-890], [Gene2], [T2], [-], [5′UTR], [protein_coding],
  [chr1:1049-1060], [Lnc1], [T3], [+], [exon], [lncRNA],
)
