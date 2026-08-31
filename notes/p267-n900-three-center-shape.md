# N900: early-shoulder thinning and concentration toward the middle

Completed source: `5f30397c5ba277fb0799fb2f7491c823de07a13d`,
`results/etop-n900-rank-width` (32M per shape pair, 800 common batches).
The `e97010a4` single-source CLI was applied unchanged: no new sampling,
family selection, or change to the primary width prediction. The constructor
was fixed on N100/N400 before reading this completed N900 archive.

The positive three-center degree-six realization exists. All 800 LOO
constructions succeed; the smallest LOO weight is 0.03031 and the largest
construction-moment error is `8.22e-15`. The empirical signed-profile negative
bin mass is 0.007805; no bins were clipped. This is an empirical positive
moment representation, not three fields or population positivity proof.

| Affine shape coordinate | N100 | N400 | N900 ± batch SE |
|---|---:|---:|---:|
| Gaussian variance fraction alpha | 0.0691064 | 0.1023710 | 0.1175059 ± 0.0183994 |
| early weight | 0.1805891 | 0.0653620 | 0.0320100 ± 0.0159030 |
| late weight | 0.3571135 | 0.4158688 | 0.3798274 ± 0.0270203 |
| relative middle gap rho | 0.2340064 | 0.2659347 | 0.3024086 ± 0.0183280 |

N900 standardized centers are `(-1.496291,-0.684928,1.186711)`; physical
centers are `(0.384733,0.481545,0.704871)`. The middle weight is
`0.588163 ± 0.023817`, continuing the point-estimate sequence
`0.462297 -> 0.518769 -> 0.588163`. The explanatory pattern is a thinning
early shoulder and increasing middle weight, **not** monotonically increasing
late weight. The earliest weight is small and uncertain; three empirical
positive atoms do not establish three population components.

In `x=N^(1/4)(p-p_ref)`, total/common-Gaussian/between-center variances are
`0.4271253 ± 0.0219792`, `0.05018975 ± 0.00579396`, and
`0.3769356 ± 0.0264533`. Relative to N400 the latter two change by
`+0.0022391 ± 0.0067922` and `-0.0435143 ± 0.0328188`; neither change is
resolved alone. The main width comparison still has no winner.

Unused standardized moments 7 and 8 are `2.7110271` and `10.8483806` versus
predictions `2.5653589` and `10.8451428`. Their residuals are
`+0.1456682 ± 0.0776455` and `+0.0032377 ± 0.1089274`;
the full-covariance score is `chi2=3.56614/2`, nominal `p=0.16812`.
This does not establish exactness of the three-Gaussian approximation.

Only the **stored** N100/N400 LOO parameter vectors from `2a824e96` were
re-expressed for comparisons; their constructions were not rerun. N900 minus
N400 has a four-coordinate affine-equality score `12.84184/4`, nominal
`p=0.012075`. Allowing alpha to vary leaves the three blur-invariant coordinates
`(early weight,late weight,rho)`, with `5.99824/3`, nominal `p=0.11170`.
Thus this recent interval does not resolve a departure from affine-plus-common
Gaussian blur. Across N100 to N900 those same gates are `121.28186/4` and
`70.38344/3`: the cumulative low-moment redistribution is not accounted for
by that blur-only orbit. These are necessary low-moment Gaussian-reference
gates, not complete functional tests; the two comparisons share N900 and
must not be combined. Consecutive increments also share N400 with covariance
`-Cov(phi400)`, explicitly saved.

Scientific card: this auxiliary changes the interpretation from width alone
to shoulder/middle redistribution in a fixed positive moment chart. It does
not identify fields, an exponent, or a fully Gaussian profile. Sources remain
three independent scale blocks with paired shapes inside each block. The next
unresolved mechanism is how the higher-moment remainder and small early-weight
limit evolve; no further collection is launched here.

Outputs: `results/p267-max-gaussian-three-center-n900/score.json` (full
800 LOO and covariance), `comparison.json` (old-LOO-only contrasts), and
`REPORT.md`. The scoring and comparison entry points are
`scripts/p267_max_gaussian_three_center.py` and
`scripts/p267_n900_three_center_comparison.py`.
