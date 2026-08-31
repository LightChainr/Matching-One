# P398: the named T4 extension repairs plus propagation, not only overlap

Parent result: `2385062d1d2e48f59d0116b56b54ffb2316c25ae`.
This is explicitly a **post-reveal deterministic mechanism exploration**.
After the preceding two-sided diagnostic identified size-four charge as a
specific omitted direction, the parent task requested this one extension.
It is not prospective validation. No other candidates, lag changes, source
changes or fitted coefficients were tried.

The only addition is
`T4=sum_{frontier blocks |C|=4} sum_{j in C} i^j` and its Kreweras partner,
on top of the previously declared triplet-incidence model. The same
stationary-L2 projection gives an **8+8-dimensional** effective dynamics.
The two exact triplet sum relations remain; all sixteen retained columns
are independent. The full named-coordinate matrix is saved.

## Actual propagation improves

| Quantity | Triplet incidence, 7+7 | Add T4, 8+8 | Full generator |
|---|---:|---:|---:|
| Crossing | .2654226141 | .2656408870 | .2656573200 |
| Minus lowest mass | 2.842890287 | 2.840410933 | 2.819658633 |
| Plus lowest mass | 1.931369756 | 1.947928395 | 1.955750138 |
| Minus t=4 relative covariance error | -4.01145% | -3.39155% | 0 |
| Plus t=4 relative covariance error | +8.06241% | +2.37679% | 0 |

The crossing's relative error shrinks from -.08835% to **-.006186%**.
The plus slow-mass error improves from -1.2466% to **-.39994%**, and its
tail overshoot drops by about 70.5%. At t=2 the plus error falls from
+2.92025% to +.81666%; the same fixed lag grid is retained throughout.
Thus this is a genuine improvement in the source's propagation, not merely
in feature/eigenvector overlap.

This is still not a uniform model-ranking claim. The earlier, smaller 5+5
pair model had plus t=4 error -2.094%, slightly smaller in absolute value
than the new +2.377%. The 8+8 model is much better in the crossing and minus
tail. No model was selected by retrospectively changing the target lag.

## The two-sided diagnostic predicted where the gain would occur

| Lowest-mode captured squared-norm fraction | 7+7 | 8+8 |
|---|---:|---:|
| Plus right | .99147931 | .99700070 |
| Plus left | .97731763 | .98016697 |
| Minus right | .99697185 | .99791898 |
| Minus left | .90981004 | .91009179 |

T4 supplies a missing plus forward geometry direction, and the improvement
survives the two-sided coupling needed to move its mass toward the true
value. The normalized plus lowest-mode residue moves from .46129229 to
.46694406, toward the full .47060221; this and the corrected mass explain
the actual long-tail improvement.

The minus adjoint mode hardly changes: about 9% of its squared norm remains outside
the geometry span. Its slow-mass error therefore improves only modestly,
from +.82392% to +.73599%. This is a concrete remaining asymmetry, not a
claim that adding more size-four labels must finish the closure.

All source moments through M^3 remain exact. The fourth moment improves in
plus but worsens slightly in minus, again illustrating why nonselfadjoint
Galerkin truncations are not a monotone variational sequence. The model is
an effective finite geometric dynamics, not an exact sixteen-state Markov
chain, a physical field count, or a Jordan construction.

## Reproduction

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /Users/lc/python-envs/research-py311/bin/python scripts/p398_width8_geometric_compression.py \
  --protocol "$PWD/analysis/p398_width8_t4_post_reveal_protocol.json" \
  --json results/p398-width8-t4-post-reveal/latest.json
```

The first invocation computed the same numbers but failed while writing
provenance because its protocol path was relative; rerunning with the
absolute protocol path saved the result. Neither model nor score changed.
No new MC, width extension, server work or further candidate search occurred.
