# Full-root one-defect gate: a single physical flip moves rank and source

## Decision

The frozen existing-fibre gate returns

```text
TWO_INDEPENDENT_DEFECT_GAIN_REJECTED
```

No new population was generated.  The calculation consumes the complete N25
axis/tilted radius-one collar fibres, imports the pooled matching root and all
componentwise Schur coefficients from the full population, and never
re-estimates a counterterm inside an edge class.

The broad Bell/rank gate contains 135,253 exact raw row classes and 8,250,462
physical pair fibres in which the same alternating site flip changes both
digital rank and canonical source Bell state.  Its first lexicographic Bell
transition has exact nonzero Schur weight `+9.6599270326e-10`.  Both sparse
kernel values happen to be zero, so this first result proves labelled Bell
motion but is not the strongest numerical-source witness.

## Stronger kernel-changing certificate

The extension was frozen before its output and additionally requires
`g16(bell0) != g16(bell1)`.  It retains 6,846 raw classes and 740,950 physical
pair fibres.  The first lexicographic nonzero transition is

```text
geometry          axis
rank              0 -> 1
source orbit       axial2
contact mask       3  (both local black arms contacted)
corner word        7
Bell               9021064 -> 2430024
g16                8 -> 0
P0                 3.2977760475e-5
P1                 4.7982147558e-5
S0                -1.5642115001e-7
S1                +1.9419147797e-9
S0+S1             -1.5447923523e-7
```

The outward rational interval for `S0+S1` is strictly negative.  Thus one
physical thermal defect simultaneously moves the topological slow variable
and the numerical canonical source.  The proposed decomposition into two
independent defects, and the automatic six-arm gain inferred from it, is
false on the exact finite state graph.

## The surviving signed functional is sharply localized

Summing every kernel-changing diagonal edge gives total signed mass

```text
-4.9488399165e-6.
```

Its source-orbit decomposition is

| source orbit | signed mass | share of total signed mass |
|---|---:|---:|
| NN other | `-3.7052834188e-6` | `74.87%` |
| diagonal | `-9.9747361552e-7` | `20.16%` |
| far | `-1.5234702796e-7` | `3.08%` |
| knight B | `-1.2215786505e-7` | `2.47%` |
| knight A | `-6.4285415345e-8` | `1.30%` |
| axial2 | `+9.2707426182e-8` | `-1.87%` |

The rank-stage decomposition is even simpler:

```text
first birth  0->1   -5.8210906659e-6   117.63% of final total
second birth 1->2   +8.7225074943e-7   -17.63% cancellation
```

Every one of the 12 rank-transition by source-orbit cells excludes zero, as
do both `S*1` row sums and all six `1^T*S` column sums.  The surviving leading
functional is therefore not a sparse accidental edge.  It is a contact-local,
first-birth-dominated four-arm channel with a smaller completion cancellation.
This is the concrete mechanism that replaces the retired independent-defect
picture.

## Consequence

The next theoretical object is no longer a generic six-arm upper bound.  It is
the signed OPE/contact counterterm carried by NN and diagonal source orbits,
followed by the separated residual after that local functional is removed.
Only the residual could still acquire a two-defect annular gain.

## Boundary

The old raw fibres retain rank, radius-one `C/B/W`, contact word and Bell
before/after, but not a canonical joint `x+y+z` component map or a single
background mask.  The nonzero transition-class sum rigorously proves that at
least one physical edge exists, but does not display one configuration-level
witness or its annular placement.  A future formal joint-map certificate can
make the witness constructive; it is not needed to preserve the mechanism
elimination above.

