# Frontier increment e: mean-clock closure versus hidden geometric transport

**2026-08-31.** This bounded continuation follows
[increment d](frontier-increment-20260831d.md). The GitHub handoff was captured
at **07:06:55 UTC / 15:06:55 Asia/Shanghai**. Two new scientific notes were
read completely, together with the P398 scientific card. No result was
recalculated. These are completed inputs to the existing
[Next Targets](../docs/NEXT-TARGETS.md), not another queue or approval sequence.

## Exact sources and integration

- **Uniform blockade and spatial variance:** `open_pr #484`,
  branch `analysis/etop-modulus-survivors-20260831`, captured head
  `d53db2f3729c722d5dd20340a8448634c75577a7`.
  [Full note](https://github.com/LightChainr/Matching-One/blob/d53db2f3729c722d5dd20340a8448634c75577a7/notes/p334-uniform-blockade-clock-semigroup.md):
  `notes/p334-uniform-blockade-clock-semigroup.md`.
  The semigroup starts at `ac9e6a8e06534f4b83bde9eefc4a7ffbf0343ea9`;
  d53db2f adds the response-variance/collision identity.
- **Current visibility versus propagation:** `branch_only`,
  `theory/p398-width8-reversible-current-control-20260831`,
  `33c6028f8722d10122b95a6c02463680cef5df5b`.
  [Full note](https://github.com/LightChainr/Matching-One/blob/33c6028f8722d10122b95a6c02463680cef5df5b/notes/p398-width8-current-source-geometry.md):
  `notes/p398-width8-current-source-geometry.md`;
  `results/p398-width8-current-source-geometry/{SCIENTIFIC_CARD.md,latest.json}`;
  definition `analysis/p398_width8_current_source_geometry.json`.

Neither citation promotes a source to main. Cross-branch discussion in #267
does not mean its scientific files have been imported into that Draft.

## 1. Uniform blockade means contain no information beyond the complete clock

For any finite monotone birth event with d selectable labels, let
`S(k)=Pr(T>k)`. Independently retain each site with probability a, leaving
blocked labels as inert dummies in the original insertion order. Averaging
over the mask and the order gives

```text
(B_a S)(k) = sum_{j=0}^k binom(k,j) a^j (1-a)^(k-j) S(j)
B_a B_b = B_(ab).
```

With a=exp(−t), the generator is
`(L S)(k)=k[S(k−1)−S(k)]`. Thus the entire mean blockade-dose response is
already closed on the full clock vector. The same conclusion holds for a
uniformly selected fixed number of blocked sites. It does not assert
closure of a smaller physical state such as (k,H2,b2).

This changes the useful discriminator. Finer mean dosage curves cannot
separate equal-clock trigger geometries. Retaining the spatial mark can.
For a uniformly selected blocked site v, write
`p_k=Pr(T=k)` and
`c_k=sum_v Pr(V_final=v | T=k)^2`. The exact singleton response obeys

```text
E_v[Delta S_v(k)] = k p_k / d
Var_v[Delta S_v(k)] = (k^2 p_k^2 / d) (c_k - 1/d).
```

The variance concerns each site's exact conditional response, not
unseparated suffix-sampling noise. In the constructed double-star versus
C4-plus-inert-site example, every mean uniform-blockade curve agrees.
At k=4, however, spatial response variances are14/625 and4/625. Conditional
final-site collisions therefore expose information lost by unmarked means.

The result is algebraic, with no new samples or network solves. Remaining-site
occupancy u is not the full-N canonical p after a nonempty prefix. If blockade
can prevent birth forever, the true mean waiting time can be infinite; finite
survival and censored-clock identities remain valid. The example is not an
observed equal-clock pair among the147 archived real prefixes.

## 2. Observing all instantaneous current does not recover propagation

In the same positive width8 process, use normalized source e and stationary
adjoint G*. Put `S=(G+G*)/2`, `J=(G−G*)/2`, and retain the unique normalized
direction r proportional to Je after removing its e component. Its maximal
initial anti-Hermitian response saturates a Cauchy bound **by construction**;
this is not a fitted optimal-observer discovery.

The existing seven-dimensional geometry already contains93.2295%/81.7318%
of the minus/plus current-source squared norm. T4 adds only.12939/.16723
percentage points. Its substantial plus-tail repair cannot primarily be
direct recovery of the initial current direction.

Indeed, the fixed pair (e,r) captures that direction completely but predicts
plus slow mass3.62846 instead of1.95575 and loses99.7429% of the t=4
correlation. The same pair fails under S. These are observable projections,
not new Markov chains.

Let P project onto (e,r) and Q=I−P; this Q is a projector, not Potts fugacity.
Because QJe=0,

```text
full initial curvature - projected initial curvature
    = ||Q S e||^2
```

for both G and S. The missing amounts are.0602551994 on minus and2.9272489452
on plus. Once the complete instantaneous current is retained, the omitted
initial feedback is reversible-force geometry. Together with the preceding
current-enabled T4 repair, this supports indirect hidden geometric transport,
without identifying one unique microscopic pathway.

All observations reuse the same generator, sources, geometry spans and fixed
lag grid. No new Monte Carlo, continuum-field count, Jordan identification
or autonomous two-state model follows.

## Current handoff, without premature completion claims

The completed objects are the uniform-blockade mean semigroup and the
current-defined observer analysis. Spatial conditional variance/site collision
and explicitly named hidden QSe transport remain meaningful distinctions;
neither requires repeating the first mean-dose or current-observer calculation.

N900 is **author-reported running, with no final score reported** at this
capture. No completion receipt or process inspection establishes more here.
The existing campaign must not be duplicated from an older “not started”
sentence.

**Delivery annotation after this capture:** the main task subsequently
completed the [P40 million-sample absolute-source Gram analysis](../results/p40-absolute-cluster/REPORT.md),
commit `a4cbf02a48c3f78ee8fb3a1e4141bd985c0bf845`. It resolves strong
global matching/source coupling but not the H4 direction difference.
It is separate from the N65/N13020k Phase-E E_top/source replay and
does not recover its missing E_top cross-products. No new source result
was imported from the two frontier branches above.

The bounded source-review subtask wrote only this note: it did not run a
scientific calculation, test or server action. The main task's separate
P40 delivery and overview/Issue updates are recorded independently.
