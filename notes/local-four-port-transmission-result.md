# A specified local pair interaction reaches original U

**Scope update:** this completed value is the **first insertion tangent**.
The subsequent [actual two-insertion witness](local-pair-two-insertion-geometry.md)
has a Q1 pole in the uncompleted tensor family. The
[regular singlet completion](local-pair-two-insertion-algebra.md) is a
different local family: its direct Q1 response is exactly zero, with
nontrivial [Q-activated interactions](regular-pair-two-site-q-susceptibility.md).
These later results do not change the first-tangent value below.

**The frozen local-interaction null is rejected on the original N25
pair.** The C4-averaged, site-average four-port insertion has

\[
\boxed{\left.\partial_\epsilon U\right|_{Q=1,\epsilon=0}
       =+0.0018155512845251097.}
\]

The exact rational enclosure of `V/A25` is

```text
lower = 194261171154185056275142086990429099 / 10000000000000000000000000000000000000000
upper = 97130585577092528137571043495214601 / 5000000000000000000000000000000000000000
```

Both endpoints are positive. They are finite-count arithmetic bounds,
not statistical confidence limits. The baseline pooled root, slope and
U/A were imported from the completed Q1 packet; no root was found again.

## The actual microscopic interaction

The [kernel](closed-source-local-four-port-pair-kernel.md) is the C4
average of the ordered unequal-pair projector `i P_[Q-2,2] i^dagger`
inserted into the vacant part of one four-edge colour vertex. The
[topology proof](local-four-port-pair-insertion.md) converts it exactly
to the original occupation mark

```text
t_x = I_(NS|EW) + (I_(NE|SW)+I_(NW|ES))/2,
S = -sum_x(t_x)/N                 at Q1.
```

Here x must be vacant, its four NN neighbours occupied, and their
components exactly the indicated two distinct pairs. The old ambient
rank of the unmodified occupation determines q/E. No virtual projector
connection is inserted into that rank calculation.

The coefficient epsilon changes each local tensor by `epsilon/N` times
this C4 average. Translation invariance makes the first-source moments
of S exactly equal to those of `-t_origin`. This is why the enumeration
uses one origin without scanning every site at each leaf. The displayed
V is in site-average units; an extensive sum-source convention would
give N times this value. No scale was fitted after the result.

This interaction is outside the original fixed-Q thermal/activity
vertex span. It is also distinct from the full seam trace: separated
winding rows support that global trace but none of these local marks,
whereas contractible local two-pair configurations support this mark
without a nontrivial global seam character.

## How it reaches U

The complete [original-observer interface](local-four-port-original-u-interface.md)
retains direct q/E source moments, covariance centering, the moving
pooled root and the thermal slope. The four prescribed terms are

| Term | Contribution to V |
|---|---:|
| centered direct thermal response | `+0.001814054032654995` |
| root movement | `+0.00027404148892317846` |
| source change of thermal slope | `+0.0000019303930697475565` |
| root change of thermal slope | `-0.0002744746301228115` |
| **sum** | **`+0.0018155512845251097`** |

The two root-related terms substantially cancel; neither was dropped.
The actual root tangents are

```text
d_epsilon h0 = -0.00042848124757398196,  h=p/(1-p),
d_epsilon p0 = -0.00007109419255347087.
```

Unlike the completed rank1-only seam packet, this local source has
nonzero raw q/E numerators. At the imported root:

| Mean | axis `(5,0)` | tilted `(4,3)` |
|---|---:|---:|
| S | `-0.0005335928112059149` | `-0.0005566767425834627` |
| qS | `+0.00044539867469842414` | `+0.0003923267884202466` |
| ES | `-0.0004568345986705169` | `-0.00040670053718612855` |

The centered q responses are positive in both geometries, so their
pooled zero moves downward. The positive original-U response is a
different question from that common root movement and is decided by
the full four-term expression above. These decompositions are not
independent evidence blocks.

## The finite birth bridge is now explicit

At a marked vacant site the activation joins two existing hypergraph
components, so `Delta cH=-1`, `Delta K=1` and

```text
Delta_x Sstar=1-ell,       ell=r(A union {x})-r(A).
```

Opposite NS|EW routing necessarily produces rank2 on activation and
has original rank0 or1, including a rank0-to-rank2 double birth at one
site. Adjacent routing occurs at every original rank. Thus there is
an explicit chain from a local colour tensor through an exterior
two-cluster routing event to a rank-transition context and the measured
original-U response. It does not require identifying the local tensor
with the already completed global seam projector.

## Scope, computation and lifecycle

- Freeze: `d7f15e68`; root score code `7c2d9fe5`, both before collection.
- Kernel proof `3c3fe12f` is integrated as `9dc3c426`; topology proof
  `ab402605` as `ba7d0632`. Both were read and accepted before GO.
- Producer `2f32d606/97f3efee` is integrated as `3b10abe1/9faa1abb`.
  Raw result `7eb1dbea` is integrated as `15372a22`.
- Each geometry had one complete 2^25 traversal for the missing local
  source crossmoments: axis1.10562s, tilted1.03976s, concurrent on the
  Mac. The full old q/E populations are retained. The single exact
  score took0.209s. No MC, cloud task, root search, source fit, support
  radius scan, old-source rescore or scientific test suite was run.
- [Raw counts and receipt](../results/p337-local-four-port-insertion/README.md)
  and [full score](../results/p337-local-four-port-insertion/score/score.json)
  preserve coefficients, input hashes, rational enclosures and the
  four response terms. This is the same exact finite N25 population
  dependency, with new prescribed crossmoments, not an independent
  stochastic confirmation.

The finite local-tensor-to-U route is complete. The C4-averaged tensor
is not a projector in one fixed cut, a free-index four-leg primary,
or a proof of continuum H4 identity. The regular one-endpoint zero
remains intact. A field claim now needs a specific long-distance or
size-dependent prediction for this fixed microscopic interaction,
not another N25 source mixture to improve its response.
