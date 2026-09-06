# No pre-existing discrete label indexes #582's remainder, and the remainder is not smooth curvature

**Date:** 2026-09-06
**Ticket:** #584 (Gate 3)
**Claim level:** C2 — reanalysis of committed productions, no new samples
**Artifact:** `results/type582-residual/latest.json`
**Scripts:** `scripts/score_type582_residual.py` (builds on `score_wasserstein_shape_flow.py`)

#584 asks whether the small, resolved, structured remainder `r_N(u)` that #582
found after its dominant direction is indexed by an already-existing discrete
label. This is the honest screen, with #584's own discipline: freeze the
dominant direction, test only pre-existing labels, and score against an exact
permutation null rather than optimism.

## What was re-frozen

`score_wasserstein_shape_flow` does not persist the per-size 9x9 delete-one
covariance `S_N`, so this module re-reads the committed histograms and re-derives
`S_N`, the affine-orthogonal residuals, a consensus `g_N` (power iteration on the
five residuals), and each transition's `r_N = v_shape` with `span{1, Q_N, g_N}`
projected out. The affine statistics reproduce #582 bit-for-bit (435 746,
192 537, 292 805, 154 876, 112 067 on 7 df). The frozen `g_N` carries 99.5%–
99.9% of each transition's `chi^2`; the amplitudes are `-2.18, -1.00, -1.71,
-0.78, -1.05` (x10^-3) on 65->130, 130->325, 85->170, 170->425, 145->290.

## The labels that were tested, and why

Every label below is computed from exact arithmetic (the parent Gaussian prime
and the committed per-size interpolation flag), never from the residuals, so the
partitions are fixed before the screen runs.

| label | classes | non-degenerate? | note |
|---|---|---|---|
| lineage | gaussian_13 / gaussian_17 / p50 | yes (3) | reference; #582 already saw it fail |
| multiplier m | 2.0 / 2.5 | yes (2) | reference |
| prime dominant-component parity | odd (3+2i, 5+2i) vs even (4+i) | yes (2) | new arithmetic label |
| prime min component | min=1 vs min=2 | yes (2) | same partition as parity here |
| prime cos4 sign | negative vs positive | yes (2) | the spin-4 sign of the parent prime |
| target interpolation flag | interpolation / extrapolation | yes (2) | N=325,425 are extrapolations |
| Smith / cyclic-noncyclic | — | **degenerate** | all eight sizes are primitive (gcd=1), the cyclic Z x Z/N quotient; noncyclic appears only at N=650/260/340 |
| deck-group / Gaussian cover word | — | **missing** | no per-size label in the tree |
| primitive homology sector | one value per orientation | too fine | 16 distinct values across 5 transitions |

The genuinely new candidate is the parent prime's **dominant-component parity**:
13 = 3+2i and 29 = 5+2i have an odd dominant component, 17 = 4+i an even one. It
coarsens lineage 3 -> 2 and was the one label that *could* have matched the
`{65,145}` vs `{85,170}` split #582 reported. It does not survive.

## The screen (covariance-identifiable, exact permutation null)

Angles are reported both covariance-weighted (primary; the common covariance is
the mean of the five transitions' `(S_base + S_target)/log(m)^2`, pseudo-inverted)
and unweighted L2 (reference, comparable to #582's note).

| label | weighted within | weighted between | weighted separation | permutation p | leave-one-out |
|---|---:|---:|---:|---:|---:|
| lineage | 13.96 deg | 17.50 deg | +3.54 deg | 0.133 | 3/5 |
| multiplier m | 15.88 | 17.40 | +1.52 | 0.300 | 2/5 |
| prime dominant parity | 16.91 | 16.71 | -0.20 | 0.700 | 2/5 |
| prime min component | 16.91 | 16.71 | -0.20 | 0.700 | 2/5 |
| prime cos4 sign | 16.53 | 16.97 | +0.43 | 0.600 | 3/5 |
| target interpolation | 15.88 | 17.40 | +1.52 | 0.300 | 2/5 |

The permutation p-value is the exact fraction of relabellings (every assignment
consistent with the class sizes; at most 30 for five transitions) whose
within-vs-between separation is at least the observed one. **No label beats the
permutation null.** `lineage` is the closest (p = 0.133) and does not clear any
honest threshold; it is the same lineage signal #582 already saw as an amplitude
effect, not an orientation effect.

## The remainder after removing `g_N` is a different object from #582's clusters

#582's `{65,145}` vs `{85,170}` clusters were read off the affine-orthogonal
shape residual `v_shape`, i.e. *before* the dominant direction is removed. The
remainder `r_N` (what #584 asks to type) is not that object. Its pairwise angles
are

```text
65->130 vs 130->325  55.2   130->325 vs 85->170  47.2
65->130 vs 85->170   40.4   130->325 vs 170->425 10.2
65->130 vs 170->425  49.2   130->325 vs 145->290 39.5
65->130 vs 145->290  27.5    85->170 vs 170->425 37.1
                            85->170 vs 145->290 54.2
                           170->425 vs 145->290 39.0
```

The tight pair is now `{130->325, 170->425}` — the two `m=2.5` transitions — not
the `{85->170, 170->425}` pair (3.9 deg in #582) that looked like a lineage
effect. So the apparent cluster structure of #582 was largely the dominant
direction's amplitude varying across lineages, and the genuine remainder has a
different, weaker, and — by the screen above — unlabelled structure. (The two
angles are different metrics: #582's are unweighted on `v_shape`, these are
unweighted on `r_N`.)

## The Taylor-curvature null: the remainder is not smooth curvature

For the two three-size lineages the second-difference curvature (unequal-step
divided difference, affine-orthogonal) is strongly resolved — statistic 1.82e6
(gaussian_13) and 1.19e6 (gaussian_17) on 7 df — but it is not a common second
direction and it does not consistently align with the remainder: the two
lineages' curvatures are 78.5 deg apart, and each curvature sits at 44.4 deg /
34.8 deg from `g_N` and at 43-80 deg from the transitions' own `r_N` (4-53%
shared variance, never dominant). Smooth second-order curvature therefore does
not provide a unified explanation of the remainder, and the remainder is not a
`d^2 Q/ds^2` artifact of a one-parameter law. `p50` has two sizes and cannot
supply a curvature.

## Verdict (decision table)

**No pre-existing discrete label beats the permutation null.** The remainder is
structured — it is not noise, and it is not smooth curvature — but its index is
unresolved on the existing five transitions. Per #584's table this is the
"downgrade the fiber" branch: keep the dominant transferable #582 direction as
the robust finite object, do not fit a free third direction or an exponent to
five residuals, and resolve the remainder with one strategically crossed
transition or a changed readout rather than with another rank.

## Not established / deferred

- **#584 step 4 (couple base/fiber to #581's typed channels) is deferred.** The
  duality-even Betti and ambient-homology tangents are exact only on the
  square-bond torus; the five #582 transitions are square-site and have no
  canonical lift, so the coupling needs a square-site decomposition that is not
  built yet. Recorded here, not done.
- A label that beats the permutation null on five transitions would still need
  #588's state-vs-memory closure before it is a state coordinate; nothing here
  reaches that point.
- The `p50` lineage's missing second-difference limits the Taylor null to two of
  three lineages.
