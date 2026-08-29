# P200 N650 context/morphism opportunity map

This post-reveal map uses only the existing 20,000-sample mixed-join batches.
It starts from the statistic actually stored by the runner,

\[
R_c=J^{\rm full}_c-\sum_{x=1}^{65}b_1(\operatorname{Inc}_{c,x}),
\]

so the isolated `C2 x C5` incidence contribution has already been subtracted
configurationwise. The raw CSV does not retain separate samplewise
`J_full/J_local` values and cannot support a different baseline fit.

## What survives the local baseline

At the frozen `p_ref`, the exact isolated-fiber means are `0.8715704195` for
black and `0.2331155417` for white. The common response across the two N650
HNF geometries is:

| typed context | local baseline / fiber | connected residual / fiber | implied full join / fiber |
|---|---:|---:|---:|
| black | 0.87157042 | -0.67716808 (SE 0.00102471) | 0.19440234 |
| white | 0.23311554 | -0.14910692 (SE 0.00054523) | 0.08400862 |

Ordinary connectivity therefore screens most of the isolated incidence
redundancy, in both typed layers. This nonzero connected remainder is a
mixed-factor interaction by construction; no chronological variable is used.

Using the exact one-orientation local fluctuation references, the scaled state
`(ES,ED,OS,OD)` is

```text
(-7.1723015, -0.0007445, -3.6803589, -0.0030854).
```

These are mechanism effect units, not sampling z scores. The covariance-aware
common-geometry test on `(ES,OS)` gives `chi2_2=945025.8310`, while the
geometry-difference `(ED,OD)` test gives `chi2_2=2.65645`, `p=0.26495`.
Equivalently, the second-minus-first geometry differences are only `0.01115`
in the even row and `0.05755` in the odd row. Thus both color-even and
color-odd static interactions survive, but no N650 embedding direction is
resolved at this precision.

## Typed H1 is an endpoint context, not a clock

The ambient-H1 state is

```text
(-1.99905, -0.00010, 0.00640, -0.00095).
```

It is global rank data and remains `O(1)`; dividing it by 65 would mix a
bounded topological rank with an extensive partition residual. Its common
subspace is almost entirely the color-even `ambient_ES` row. The matching-odd
common row has `z=0.875`, and ambient `(ED,OD)` gives `chi2_2=0.52444`,
`p=0.76934`. This is a convention-labelled endpoint defect under the frozen
representative-displacement lift. It is correlated with the primary stream
and is not a second evidence block.

## Why this is not path/state memory

`ED/OD` compare two static HNF embeddings evaluated with shared counters. They
do not compare `R2 then R5` with `R5 then R2`; the endpoint joins commute.
Moreover, the 20k schema stores neither the four intermediate states
`h0,h2,h5,h25` nor an activation time, marked lineage, or filtration
increment. Consequently:

- nonzero `ES/OS` identifies a static connected mixed-factor response;
- null `ED/OD` supports a geometry-common response at N650;
- neither result identifies or rejects chronological memory.

Calling the large `OS` signal path memory would therefore be a type error: it
is matching-odd under the typed color layer, but it has no path index.

## Conditional morphism freeze

The two-geometry agreement permits a post-reveal parameter freeze for a new
N650 HNF embedding with the same typed observable:

```text
E = -53.707875,  O = -34.323975.
```

For a genuine scale challenge, freeze the connected residual densities

```text
rho_B = -0.6771680769 +/- 0.0010247146 per source fiber,
rho_W = -0.1491069231 +/- 0.0005452309 per source fiber.
```

An explicit unrun N1300 pair keeps relative factors `1+i` and `2-i`:

| source | N260 | N650 | N1300 final |
|---|---|---|---|
| `11+3i` | `8+14i` | `25-5i` | `30+20i` |
| `7+9i` | `-2+16i` | `23+11i` | `12+34i` |

Under the unverified hypothesis that the residual is extensive in the 130
source fibers and geometry differences remain zero, the frozen state is

```text
(ES,ED,OS,OD) = (-107.41575, 0, -68.64795, 0),
R_B = -88.03185, R_W = -19.38390.
```

The fit-only SEs on the last two color residuals are `0.13321` and `0.07088`.
No scale discrepancy can be estimated from one N, so these numbers are a
future falsifiable morphism-density challenge, not an established scaling law
and not authorization for new production.

## Opportunity surface

| axis | current result | scientific status |
|---|---|---|
| local incidence | removed exactly in the stored residual | exact observable semantics |
| typed color | black and white residual densities both nonzero | mixed-factor interaction |
| HNF embedding | ED/OD compatible with zero | geometry-common opportunity |
| ambient H1 | common-even endpoint defect near -2 | convention-labelled secondary |
| join chronology | absent from raw schema | path/state memory not identifiable |
| scale | one N only | conditional N1300 freeze, unvalidated |

The five-line scientific card is stored beside the machine-readable result.
No new acquisition was started.
