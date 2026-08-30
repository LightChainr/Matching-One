# Complete F3 constraint curves and the N65 geometry crossing

Status: zero-new-sample analysis of the P334 N65 20k projective-birth archive.
The exact output is a curve engine and covariance-sufficient Bernstein archive;
the root result is exploratory design evidence.

## The sparse event tuple already contains the whole curve

For each path, occupation count `k` and nonzero `alpha in F3^2`,

```text
k < tau1:             T_alpha(k)=1,
tau1 <= k < tau2:     T_alpha(k)=1[alpha dot ell = 0],
k >= tau2:            T_alpha(k)=0.
```

The zero twist is exceptional and stays `T_00(k)=1`.  Averaging paths gives
the complete microcanonical coefficients

```text
T_alpha(p)=sum_k c_alpha(k) B_(N,k)(p)
           =P0(p)+L_kernel(alpha)(p),   alpha != 0.
```

Derivatives are analytic without finite differences:

```text
T_alpha_prime(p)
 =N sum_(k=0)^(N-1) [c_alpha(k+1)-c_alpha(k)] B_(N-1,k)(p).
```

All eight nonzero twists are retained, although proportional `alpha` and
`2 alpha` have the same kernel and therefore give four identifiable nonzero
constraint curves.  Every nonzero curve is monotone decreasing; `T_00=1` has
zero derivative.  Consequently a zero/nonzero crossing exists only at the
trivial endpoint.  The informative internal crossings are differences between
nonzero sectors or their exact projective characters.

The machine artifact stores the aligned per-batch degree-65 coefficients for
`P0`, `P2` and all four F3 line bins.  Those are sufficient statistics for
arbitrary-p, derivative and cross-p covariance.  A 32-point central grid also
contains the full within-p 48x48 covariance across both orientations, nine
twists, nine derivatives, three characters and three character derivatives.

## Exact complement and D4 transport

Alexander/complement reversal acts as

```text
P0_dual(p)       = P2(1-p),
P2_dual(p)       = P0(1-p),
L_line_dual(p)   = L_line(1-p),
C_dual(p)        = C(1-p),
C_dual_prime(p)  = -C_prime(1-p).
```

Thus a nonzero twist transports to

```text
T_alpha_dual(p)=P2(1-p)+L_kernel(alpha)(1-p),
```

not to `1-T_alpha`.  This distinction prevents the common rank-one piece from
being assigned the wrong complement parity.

In the normalized character basis `[H,A,D]`, exact projective D4 acts by

```text
S quarter-turn:  (H,A,D) -> ( H,-A,-D),
x reflection:    (H,A,D) -> ( H, A,-D).
```

Microcanonical partition, zero/proportional-twist, monotonicity, D4 and
complement gates all pass; the largest numerical residual is `2.50e-15`.

## A parameter-free physical-geometry crossing

In the central window `[0.45,0.75]`, roots were computed for each character on
each physical N65 Gaussian orientation and for the second-minus-first curve.
Only one candidate survives the same bracket in all 20 leave-one-batch
reconstructions:

```text
H_(7,4)(p_cross) = H_(8,1)(p_cross),
p_cross = 0.573633326,
d/dp [H_(7,4)-H_(8,1)] = 0.0602462.
```

The F3 H weights and the equality contain no fitted amplitude.  At the root
both curves equal `0.1471914013`.  Exact complement transport predicts the
dual matching-lattice partner

```text
p_cross_dual = 1-p_cross = 0.426366674,
slope_dual = -0.0602462.
```

The sign pattern on the reused block is

| p | H second-minus-first | SE |
|---:|---:|---:|
| 0.500000 | -0.00341247 | 0.00111025 |
| 0.550000 | -0.00140576 | 0.00122653 |
| 0.573633 | 0 by construction | 0.00137529 |
| 0.592746 | +0.00110298 | 0.00150621 |
| 0.650000 | +0.00285719 | 0.00153531 |

The root jackknife SE is `0.02281`; its distance from the project reference is
only `0.01911`.  Therefore the 20k archive does not resolve a physical root
split from `p_ref`, even though it supplies a stable geometry-selector bracket.
Under simple inverse-sample scaling the present variance implies approximately
114k samples/shape for a 2-sigma separation and 256k for 3 sigma.  These are
power numbers, not a production request; no new sampling was performed.

The other aggregate central roots have leave-one-batch bracket survival from
55% to 85%, or disappear entirely, so they remain unselected.  The 20/20 H
crossing is the only defensible frozen target carried forward.

## Scientific use

This result supplies a genuine p-dependent explicit twist observable without
adding a synthetic source to the Monte Carlo engine.  It changes the next
experiment from “measure diagonal odd again at one p” to “score the frozen H
geometry crossing and its complement partner on independent counters or a
child geometry.”  It does not yet establish a continuum defect, a universal
crossing value or an independent charged response.

Reproduce with:

```bash
python3 scripts/score_f3_twist_curves.py \
  --births results/local-20260830/P334-projective-birth-N65-smoke/n65_20k.births.csv \
  --metadata results/local-20260830/P334-projective-birth-N65-smoke/n65_20k.metadata.json \
  --json results/local-20260830/P337-F3-twist-curves/score.json \
  --markdown results/local-20260830/P337-F3-twist-curves/REPORT.md
```
