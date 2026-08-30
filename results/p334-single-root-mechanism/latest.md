# Why each projective orbit has one simple balance root

Status: `exact_single_root_criterion_and_bounded_ulc_audit`.

## Exact theorem

```text
d_k=(k+1)A_(k+1)-(N-k)A_k=(N-k) binom(N,k) (q_(k+1)-q_k)
```

After deleting zeros, q_(k+1)-q_k has one sign change from positive to negative, with both signs present.
Then A(p) has exactly one simple critical point in 0<p<1.

With t=p/(1-p), factor A'(p)=(1-p)^(N-1) P(t). The positive coefficients of P have strictly lower degrees than all negative coefficients. Writing P=P_plus-P_minus, P_minus/P_plus is strictly increasing because its log-derivative is the difference of two degree averages with disjoint ordered supports. It crosses one exactly once and transversely.

Thus strict unimodality is stronger than necessary: a two-layer plateau merely inserts a zero derivative coefficient and does not destroy the unique simple root.

## Exact audit

The focused six-HNF atlas plus N13/N17 contains 16 orbit rows; all are strict single peaks with strict ULC/ratio decrease: `True`.

The broader scan covers 83 honest connected HNFs, 240 fixed-line sequences and 217 stabilizer-orbit sequences.

| property | fixed lines | line orbits |
|---|---:|---:|
| contiguous support | 240/240 | 217/217 |
| weak single peak | 240/240 | 217/217 |
| strict single peak | 228/240 | 205/217 |
| strict ULC on support | 240/240 | 217/217 |
| strictly decreasing ratios | 240/240 | 217/217 |

The first failure of *strict* single-peakedness is already N=6 at `[[2, 0], [0, 3]]`, group `[[1, 0]]`: q on support is `['1/5', '3/5', '3/5']` with modes `[3, 4]`. It remains strictly log-concave and its ratio sequence remains strictly decreasing.

No log-concavity, ratio-monotonicity or one-sign-change counterexample occurs through N12. Therefore the exact theorem certifies a unique simple root for all 217 orbit rows without numerically solving their polynomials.

## Structural boundary

For fixed primitive line ell, 1{H1=ell}=1{ell subset H1}-1{rank(H1)=2}. Both terms are increasing events; rank two is contained in the line-containing event.

The fixed-line rank-one family is order-convex in the Boolean lattice because homology images grow under site inclusion.

But `[[2, 0], [0, 2]]` at N=4 already has a rank jump 0->2 when site 0 is added to sites [1, 2]. Therefore the literal site-matroid route fails: one-site rank increments need not be at most one.

A difference of nested increasing events need not have a log-concave layer density. On four elements, U={sets containing 1} union {sets containing {2,3,4}} and V={sets containing 1 with size>=2} give U\V={{1},{2,3,4}}, whose normalized layer sequence has a zero valley and is not log-concave.

## Theorem/conjecture ladder

- **THEOREM**: one sign change of normalized layer differences implies one simple root
- **EXACT_BOUNDED_EVIDENCE**: All 240 fixed-line and 217 stabilizer-orbit sequences in the honest connected HNF N<=12 scan are contiguous, strictly ULC on positive support, and have strictly decreasing adjacent ratios.
- **CONJECTURE**: Projective rank-one ULC conjecture: every honest finite square-torus quotient and fixed primitive line has a contiguous ULC normalized layer sequence, hence a unique simple balance point.
- **PROOF_TARGET**: Show the order-convex line stratum has a normalized matching/Lorentzian rank-generating property; ordinary matroid and nested-upset log-concavity do not suffice.

## Boundary

- The unique-root implication is proved; projective-line ULC for all quotients is conjectural.
- Strict single-peakedness itself already fails at N=6 through a two-layer plateau, but the one-sign-change criterion survives.
- The matroid obstruction rules out the literal site-ground rank model, not every possible auxiliary lift.
- No Monte Carlo sample, Huawei production, new PR, or merge is used.
