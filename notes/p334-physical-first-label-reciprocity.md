# P334: real first-label information recovered from the existing final-site clock

## Result

The two already solved N425 prefixes both have **zero original direct sites**,
so their binary direct/safe first-step Doob floor is zero. Observing the actual
first insertion label nonetheless removes **3.12935% / 3.68332%** of their
conditional clock variance. The whole temporal response is strongly, but not
exactly, two-mode: B assigns 2.58 times A's trace fraction to its second mode.

| Exact conditional scalar / numerical Gram readout | A: 43042514269 | B: 43042505280 |
| --- | ---: | ---: |
| Remaining labels d | 173 | 173 |
| Original direct count h | 0 | 0 |
| E[T] | 17.73237780 | 20.77877866 |
| Var(T) | 97.52952154 | 147.04286285 |
| Var(E[T given first label]) | 3.052035868 | 5.416052259 |
| Fraction of Var(T) explained by first label | 3.12934568% | 3.68331530% |
| Survival-Gram first-mode trace fraction | 98.53652573% | 96.28300016% |
| Survival-Gram second-mode trace fraction | 1.41094737% | 3.64024255% |
| First two modes combined | 99.94747310% | 99.92324271% |

The scalar quantities are computed with exact rational arithmetic; the table
shows decimals. The Gram spectrum uses double precision, with the complete
matrix and curves saved. Its small eigenvalues are **not** exact-rank claims.
This is one thin readout of the old selected-prefix law, not a new ensemble or
an extension of the 40-batch production result to all prefixes.

## Coefficients and one recurrence

The source is `1c06230b`, specifically `true_safe_counts` in
`results/p334-contracted-full-clock/full_physical_birth_clock.json` and the
173 per-site `pivotal_count_by_prior_size` arrays in
`results/p334-exact-marked-birth/marked_birth_{counter}.json`.
The source script defines their entries as integer forced-off minus forced-on
safe counts over the other 172 labels. Thus they are exactly b_v(k-1) in the
first/last-label reciprocity theorem at `31c17d48`:

\[
c_v(0)=0,\quad c_v(k)=I_{k-1}-c_v(k-1)-b_v(k-1),\qquad
s_v(k)=\frac{c_v(k)}{\binom{d-1}{k-1}}.
\]

Here s_v(0)=1, and

\[
m_v=\sum_{k=0}^{d}s_v(k)=1+E[T_{\rm child}(v)],\qquad
\Gamma_{kl}=\frac1d\sum_v(s_v(k)-s(k))(s_v(l)-s(l)).
\]

The exact rational scalar innovation is
`mean_v (m_v-E[T])^2`; the total clock variance comes from
`E[T^2]=sum_k(2k+1)s(k)`. The integer recurrence recovers
`sum_v c_v(k)=k I_k` at every horizon. No physical event, child continuation,
partition DP or permutation is solved again.

## What the first label reveals geometrically

In A, the eight tied port-0 sites 6,27,140,162,251,274,296,409 reduce the
conditional mean to 11.69773835 when inserted first. In B, site121 reduces
it to 6.59390871; sites8,166,277 each give 7.18272346. These are conditional
means **after choosing the first label**, not the previously reported times
conditioned on which label wins at the end.

All exactly irrelevant final-winner sites remain informative as first labels:
the 90 in A give mean18.63046758; the 46 in B give21.65936039. For any irrelevant
label, `m_inert=1+d E[T]/(d+1)`: choosing it first spends one step without
advancing the event. This class contributes 13.7482% / 3.80687% of the total
first-label innovation. Even their shared clock-delaying effect is not zero.

The four existing physical roles (port0, port1, interior, outside core)
explain only 35.2400% / 21.8961% of the full label innovation. Thus a role-only
observer discards substantial information within the named network roles,
especially in B. No new site classification was fitted.

The two numerical Gram spectra begin

```
A: 0.1037989908, 0.001486300759, 0.00005472998950, ...
B: 0.1564448131, 0.005914824679, 0.0001184700022, ...
```

The real examples retain the toy lesson: a common zero binary floor does not
specify full-first-label information, and the allocation among temporal
response directions differs. The examples do **not** have equal unmarked
clocks, so they are not a second physical isoclock counterexample.

## Dependency and scope

The reciprocity transform is invertible in time. Consequently full first-label
curves and the full time-resolved final-winner law are equivalent information;
their response ranks agree algebraically. The new Doob/Gram interpretation
does not add an independent evidence block to `1c06230b`. Marginal winner
probabilities or their collision alone would be insufficient.

Both prefixes are the same second-orientation N425, k0=252, seed20260831430425,
age10, ell=(12,-19) selected physical states used previously. First-label
innovation is conditional variance reduction, not intrinsic temporal memory,
an H4 paired bound, continuum-state count, or a population scaling claim.

Run (single-thread local, about0.24s):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
/Users/lc/python-envs/research-py311/bin/python scripts/p334_physical_first_label_reciprocity.py
```

`score.json` retains exact per-label means/innovation fractions and immutable
source hashes. The two NPZ files retain the ordered labels, all174 horizons,
conditional survival curves, numerical Gram matrices and eigenvalues.
