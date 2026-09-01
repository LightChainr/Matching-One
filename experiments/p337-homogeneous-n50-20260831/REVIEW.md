# P337 homogeneous N50 independent scoring review

**PASS.** The fixed `(5,5)/(1,7)` homogeneous N50 target reproduces with an independent direct-Bernstein calculation. No table or frozen scorer was modified; no new point, model, population, or source was evaluated.

## Method and integrity

- `CONTRACT.md` and `score.py` are byte-identical to freeze `10c666b65566b25ddb8eaa02219947a9c5a261f2`; the pinned interval-backend hash also agrees.
- Both complete integer tables have exactly `2^50 = 1125899906842624` configurations. Every K layer sums to `binomial(50,K)`; rows are unique `(K,q)` entries with integer `count,sum_S`.
- Input SHA256: first `a6f2ae1a31eca7c8e08e271bf28639f68369a3585ef8384d04d776bc78d933e7`; second `d0f7bb183a5bfbdcb00979b22e92a59775ffb70d96b8859eec0be5d2e9225df5`. These are exactly the inputs identified by the frozen primary result.
- The independent script [review.py](review.py) (repository-portable copy of the executed temporary script) directly differentiates `p^K(1-p)^(50-K)` through second order for `1,q,E,S,qS,ES`. It independently normalizes each geometry and uses `J_f=<fS>-<f><S>` and its p derivative. It imports neither the primary scorer nor its interval backend and does not use the logit K-moment formulas.
- Root and denominator motion are included. The four p-coordinate terms are converted to the frozen z coordinate by the chain rule; the two root terms change separately while their sum is invariant.

## Numerical reproduction

| Quantity | Independent Decimal result |
|---|---:|
| Common root p | 0.5927594013067592658715481382408770 |
| U | 1.0615603876876550536637165689584383 |
| V | +0.05434578266955829908314715506295295 |

The maximum difference of the reported core values at Decimal precision 120 versus 160 is `3.16073e-116`. All 12 checked values lie inside the primary exact outward-rounded Fraction intervals: p, h, slope, U/A, V/A, root source motion, all four response terms, and both geometry S means.

In particular, the primary rational lower bound is

`V/A >= 1885252764556849639086318759550302603 / 10^40 > 0`,

with serialized interval width `10^-40`. The rational interval establishes the sign. The Decimal calculation is a separate numerical cross-check, not an additional strict interval certificate or independent statistical evidence.

## Decision and limits

The frozen positive-sign continuation prediction is **not rejected**. The finite-point null `V=0` is excluded at this specified N50 law and fixed source. This does not confirm a mechanism or prove continuation throughout epsilon, other sizes, or an asymptotic regime. This audit checks the supplied sufficient-statistic tables and scoring mathematics; it does not constitute a separate enumeration proof of the producer.

Machine-readable reproduction: [results/review.json](results/review.json). Primary result SHA256: `e4af39f02bfaf4a7ac0eb41af723dcf9aa7498d64009eecb167d1ab3bc6515ba`.
