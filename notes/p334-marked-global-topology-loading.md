# The first birth changes the sign of the marked topology source contrast

On the same accepted global pairs, adding the first-birth term reverses both
integrated source contrasts at each size. A source's positive completion-F2
loading therefore does not determine its contribution to complete A_top.
The newly evaluated quantity is the physical final-source-marked full
topology response, not a renamed F2 proxy.

These are point estimates from the existing fixed populations; source signs
and size changes have not been independently established. The original
20-batch vectors are handed to the single global covariance coordinator.

## Exact source readout

For each accepted R1 prefix, let s be either the original H2 direct set or
its collective complement, with completion law p_s(j), winning probability
`pi_s=sum p_s(j)` and source-weighted delay `tau_s=sum j*p_s(j)`.
The already-fixed full observable is `A_top=g(K1)+g(K2)-1`, where
`g(k)=P(Binom(N,p)>=k)`. Since K1 is known from the ordered prefix,

\[
 A_s(p|\mathcal F)=E[A_{\rm top}(p)1_{\text{final source}=s}|\mathcal F]
 =\sum_jp_s(j)g(k_0+j)-\pi_s[1-g(K_1)].
\]

The first term is a nonnegative orientation-level completion contribution.
The second is the first-birth subtraction; its orientation contrast can
have either sign. Both are stored on the **same** accepted pairs, so their
comparison is not caused by changing the gate between terms.

For the integral,

\[
 A_s^{\rm int}=\frac{(d+1)\pi_s-\tau_s-K_1\pi_s}{N+1}.
\]

The direct law `p_D(j)=h*S(j-1)/(d-j+1)` implies
`(d+1)*pi_D-tau_D=h*mu`, where `mu=sum_{k=0}^{d-1}S(k)=E[T|F]`.
Consequently

\[
 \boxed{A_D^{\rm int}=\frac{h\mu-K_1\pi_D}{N+1}}.
\]

This exposes the relevant coupling: the first-birth time is weighted by
the future source's winning probability. It cannot generally be replaced
by a product of separately averaged K1 and pi_D. No site knockout or causal
intervention is performed.

## Same-gate integrated mechanism

All entries are first-minus-second divided by the original delta_cos4.

| Size/source | Completion F2 | First-birth subtraction | Complete marked A |
|---|---:|---:|---:|
| N325 direct | +0.0005547302 | +0.0006169392 | **-0.0000622091** |
| N325 collective | -0.0001240242 | -0.0002681859 | **+0.0001441617** |
| N425 direct | -0.0003571092 | -0.0005986182 | **+0.0002415090** |
| N425 collective | +0.0008358653 | +0.0012763700 | **-0.0004405047** |

Both marked source directions reverse relative to the same-gate completion
contrast at both sizes. This is the precise reason a completion-only account
can assign the wrong source direction to the full topology integral. The
complete integral reads `1-(K1+K2)/(N+1)`, so it includes the first birth's
alignment with the final source, not only the second birth's delay.

At the canonical reference p, the corresponding full marked direct/collective
contrasts are `(+0.0011283058,+0.0005212270)` at N325 and
`(+0.0001583048,+0.0006799754)` at N425. Their competition is readout-dependent;
the original canonical and integrated vectors stay jointly paired.

## Return to the complete safe-global hybrid

The global gate is exactly the sufficient policy in
`notes/p334-global-two-birth-loading-policy.md`: both orientations have rank
at least one, and every R1 clock belongs to the original `exact_pair`
acceptance. Both R2 are identity replacements. Any R0 or original clock
failure retains the entire raw global paired vector.

Only globally accepted R1 directions enter the marked direct/collective
columns. All other full A values enter a signed remainder, including known
past R2 contributions and entire fallback pairs. There is no partial source
classification of a fallback and no missing-stratum zero substitution.

| Size/readout | Marked direct A | Marked collective A | Remainder A | Complete hybrid A |
|---|---:|---:|---:|---:|
| N325 canonical | +0.0011283058 | +0.0005212270 | -0.0017692713 | -0.0001197385 |
| N325 integrated | -0.0000622091 | +0.0001441617 | +0.0000267190 | +0.0001086717 |
| N425 canonical | +0.0001583048 | +0.0006799754 | +0.0092690377 | +0.0101073179 |
| N425 integrated | +0.0002415090 | -0.0004405047 | +0.0011420301 | +0.0009430345 |

The accepted-R1 source sum does not dominate the full contrast automatically.
N325's canonical source sum is largely opposed by the remainder; N425's
canonical point contrast is mostly in that remainder. The latter is not
identified as an R0-only mechanism: it also contains R2 and clock-policy
fallback, to be separated by the global nine-state readout.

N325 has 9055 accepted pairs involving R1, 1702 both-R2 identities,
9207 pairs containing R0 and 36 clock-policy fallbacks. N425 has
8903, 1572, 9413 and 112 respectively. Marked orientation counts are
6512/6529 at N325 and 6412/6439 at N425. All 20000 original paired counters
remain in each denominator.

## Scientific card and common covariance handoff

- Changes: the original completion mark is now attached to full A_top. Its
  first-birth coupling reverses the same-gate integrated source point
  contrasts, demonstrating a concrete distinction from the old F2 loading.
- Observer/sector: canonical full topology and its p-integral, marked by
  the final source only on globally accepted R1 directions; signed full
  remainder preserves the rest of the observable.
- Source: complete original births `9c495ab1`, the same original40k counters
  as exact-clock archive `0d1e586d`. No new MC, DP, network or source refit.
- Dependency: each size retains original batches0..19 with1000 paired
  counters. All means, source terms and first-birth subtractions are one
  dependency group with the global covariance analysis; no source p-values
  or independently combined errors are supplied here.
- Lifecycle: fixed full-birth replay -> saved conditional clocks -> global
  paired gate -> exact marked conditional A -> common-batch covariance.
- Next readout: jointly assess first-birth versus completion contributions
  and the nine-state remainder. A physical source label on an accepted R1
  slice is not an explanation of the entire global signal or a field count.

`results/p334-marked-global-topology-loading/score.json` contains source
hashes and, for each size, the 20-by-12 matrix
`joint_20_batch_means_orientation_readout_source`. Its labels are
orientation(first,second) x readout(canonical,integrated) x
source(original_H2_direct_A,collective_A,remainder_A).
`safe_global_hybrid_A_20_batch_means` has order
`first_canonical,first_integrated,second_canonical,second_integrated`.
The aligned positive-F2 and first-birth-subtraction 20-by-8 matrices use
the same ordering with just the first two source labels. Source-clock
vectors preserve pi_D,pi_G,tau_D,tau_G per orientation.

Component additivity differs from the direct hybrid computation by at most
`7.78e-16` per record. Only stored coefficients were read, and no old
reliability calculation or stochastic production was rerun.

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/p334_marked_global_topology_loading.py
```
