# Independent review: specified oracle-centred iid estimator

**PASS — no required correction.** The variance formula, locked input law,
original-U conversion and integer sample thresholds agree. This review did
not import or run `score.py`, change `result.json`, refine a root, enumerate
configurations or draw samples.

## Formula and target

In each geometry at its law's known pooled root, define
`mu_g=E_g K`, `I1=1{rank=1}` and `X_g=(K-mu_g)I1`.
Then `c_g=E_g X_g=Cov_g(K,I1)`. For independent iid samples with **n in
each geometry**, the difference of sample means is unbiased for
`theta=c_axis-c_tilted`, with

```text
Var(theta_hat) = (Var(X_axis)+Var(X_tilted))/n,
n_for_SNR_s = ceil[s²*(Var(X_axis)+Var(X_tilted))/theta²].
```

The denominator is theta², not the individual rank-one probability.
The requested equal-per-geometry allocation is retained; n is not the
combined two-geometry count.

Because `E=q²=1-I1`, `partial_log(h) E=-Cov(K,I1)` and
`partial_log(h) q=Cov(K,q)`. Therefore the original ratio is

```text
U/A25 = -theta / [DeltaCos4 * mean_g Cov_g(K,q)],
DeltaCos4=1152/625.
```

The same h-to-p Jacobian cancels from numerator and denominator. The minus
sign is essential: star theta is positive and U is negative; drop reverses
both signs. Each geometry retains its own partition normalizer.

## Independent calculation

`verify_independent.py` reads the 269 axis and 247 tilted `(K,g,q,count)`
rows. Counts, positive multiplicities, unique row keys and every
`sum_{g,q} count=binom(25,K)` are checked. The two CSVs, original m64 root
result and arithmetic vendor are byte-identical to fixed input commit
`cae9c8997b5994c218bfe060f75656137f745755`. The target code hash also matches
the hash recorded in `result.json`.

Using independent Decimal arithmetic at **120 and 160 significant digits**,
the verifier evaluates only the saved root lower endpoint, midpoint and
upper endpoint. It computes directly, row by row,

```text
E[X²] = sum_rank1 w(K)*(K-mu)² / Z,
Var(X) = sum_all w(K)*((K-mu)I1-c)² / Z,
Var(K | rank1) = sum_rank1 w(K)*(K-E[K|rank1])² / sum_rank1 w(K).
```

It does not copy the scorer's second-moment expansion. All cell quantities,
theta, the variance sum, both SNR budgets, the original q/E derivatives,
pooled denominator and original U ratio pass **420 saved rational-enclosure
checks**. The evaluated root endpoints have the required opposite pooled-q
signs. Decimal values are converted to exact Fractions for containment
comparisons; displayed floating-point midpoints are never used for a check.

Selected independently reproduced values:

| Law / geometry | E[X²] | Var(X) | E[K|rank1]-E[K] |
|---|---:|---:|---:|
| star / axis | 3.476363790877539e-14 | 3.476363790877539e-14 | 0.0001513826343874195 |
| star / tilted | 2.637719015161550e-19 | 2.637719015161550e-19 | 1.194721403576091 |
| drop / axis | 1.955114322320590e-13 | 1.955114322320583e-13 | -6.348479801408117 |
| drop / tilted | 2.758552988032879e-19 | 2.758552988032879e-19 | -1.971627537148633 |

The star axis has appreciable conditional K variance but almost cancelling
conditional mean displacement; the selected estimator pays for that
within-rank variance. This explains why a first-rank-one-occurrence budget
does not measure this estimator's signal-to-noise requirement.

| Law | theta | Original U/A25 | Exact minimum integer n per geometry for SNR 3 |
|---|---:|---:|---:|
| star | 1.435641269315531e-19 | -6.232603901919505e-21 | **15180258044365917963223558** |
| drop | -2.639984219178268e-14 | 1.146027183860387e-15 | **2524716496431136** |

For SNR 1 the corresponding integers are `1686695338262879773691507`
and `280524055159016`. In all four cases, exact Fraction ceilings of the
saved **lower and upper** budget bounds coincide; the Decimal-derived
ceilings at all three root points and both precisions agree. Thus the
integer result does not depend on a rounded double or an unresolved
fractional crossing of an integer. The original rational interval arithmetic
supplies the certificate; this Decimal calculation independently cross-checks
the defining row-square quantities rather than claiming a second formal
interval proof.

## Scope and reproducibility

This is a second-moment budget for the stated iid estimator with oracle
root, separate exact means K and exact denominator. It is not a 95% coverage
or power guarantee, an all-estimator/all-algorithm lower bound, or a cost
for optimized allocation, conditional sampling, importance sampling or
other control variates. Unknown means/root/denominator and autocorrelation
are excluded by the definition, not shown to be free in practice.

`review.json` records both precisions, all three reused root points, cell
readouts, interval-check names, provenance hashes and the single verification
run (about 0.087 seconds). Reproduce with:

```bash
/Users/lc/python-envs/research-py311/bin/python verify_independent.py
```

Original scorer SHA256:
`0b6084d5119af972084d0626a195fe85fa101de82c022eb8fd906cd48e6e3387`.
The result and input SHA256 values are preserved in `review.json`.
