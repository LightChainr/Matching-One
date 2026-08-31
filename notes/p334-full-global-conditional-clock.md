# P334: complete topology cancels the R1 membership step

## Scientific result

The large membership-step variance of the previously gated R1 second-birth
observable does not survive in complete `A_top=F1+F2-1`. On the same 20,000
paired paths per size, the fraction of integral variance explained by the four
old R1 flag pairs changes from **99.849%/99.864%** to
**0.0185%/0.00511%** for the baseline full observable (N325/N425). For the
safe conditional hybrid it is 0.0223%/0.00545%. These small empirical fractions
have substantial relative batch uncertainty; they establish neither an exact
zero nor a new variance law.

The full observable also retains a useful, differently sized conditional-noise
component. The paired residual second moment divided by baseline individual
variance is **2.317% +/- 0.0392pp / 2.323% +/- 0.0480pp** at `p_ref`, and
**4.632% +/- 0.1048pp / 4.611% +/- 0.1223pp** after integration. Observed
safe/baseline individual variance ratios are respectively
`0.97655/0.97765` and `0.95405/0.95539`. Residual-noise fractions and empirical
variance reductions are separately reported, not equated exactly on a finite
archive; neither includes implementation cost or implies an end-to-end speedup.

## Observer, source and replacement rule

The sources are final full births at
`9c495ab13e65f2bc93dc0849ee3b73f88724c4b1` and final conditional clocks at
`0d1e586dafbade5e7d1f9bfc598170d0c881e337`. Both consume the original
`e81dd59` paired permutations, not new Monte Carlo. The join is by
N/batch/counter/orientation. N325 and N425 are separate blocks; within each size
all observers retain their common 20 batches of 1,000 paths. This is one source
dependency group per size, shared with the earlier R1/source readouts.

For `g(K,p)=Pr[Bin(N,p)>=K]`, use

```
A(p) = g(K1,p)+g(K2,p)-1
integral A = 1-(K1+K2)/(N+1).
```

Conditional replacement is allowed only when **both checkpoint ranks are at
least one and every R1 clock is exact**. Then K1 is measurable from the common
ordered prefix, R2 needs no replacement, and each R1 F2 becomes its conditional
mean. Any R0 or solver failure keeps the entire original pair. Accepted pairs
with an R1 replacement number 9,055/8,903; R0 blocks 9,207/9,413; solver
fallbacks block 36/112; both-R2 pairs number 1,702/1,572. No partially solved
fallback clock is used.

The safe-minus-baseline paired mean shifts are:

| Readout | N325 | N425 |
|---|---:|---:|
| canonical | -0.000288845 +/- 0.00110700 | 0.000847299 +/- 0.00102472 |
| integral | -0.000002222 +/- 0.000166747 | 0.000180498 +/- 0.000142793 |

Mean uncertainty uses the original joint batch covariance. No stratum errors
are independently summed, and no high-dimensional covariance is inverted.

## Exact origin cancellation, measured in the full archive

Fix `alpha=1-2*p_ref=-0.18549210158`, without tuning. In the integral
H4-normalized contribution of rank r, subtract
`alpha*(1[R_first=r]-1[R_second=r])/delta_cos4`. The three subtractions sum
to zero path by path. The safe integral stratum variance sum changes from
`0.0869747 -> 0.0122105` at N325 and `0.0622048 -> 0.00753563` at N425;
the total observable and its variances `0.0146283/0.00906587` stay exactly
unchanged. Covariances supply the difference between the sum of stratum
variances and the variance of their sum. This makes the old large gated
constant step an observer-origin effect rather than 99.85% removable global
noise.

## Reproduce and next physical readout

```
python3 scripts/p334_full_global_conditional_clock.py --full-commit 9c495ab13e65f2bc93dc0849ee3b73f88724c4b1
```

The score saves source hashes, baseline/safe and R0/R1/R2 joint means, full
20-batch LOO covariance, individual covariance for the origin transformation,
and the precise safety counts. This scorer performs zero DP, zero birth
replays and zero new MC. The next available readout is the correctly marked
full-A direct/collective/remainder split (including the first-birth debt), not
the old F2-only source split. A first thermal moment can then expose the W²
clock component absent from the complete zeroth integral. These will be joined
on these same batches, not treated as independent confirmations.
