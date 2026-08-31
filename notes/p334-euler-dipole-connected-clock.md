# Thermal dipole, center displacement and birth-clock spread

The N425 source-minus/geometry-difference dipole is weak because a resolved
center-displacement contribution is partly canceled by unequal marginal
broadening. A distinct, stronger rank-covariance response remains after
subtracting shared order-statistic timing.

This analysis uses the original common-label fork block. It combines the final
thermal moment readout `9059776d` with the full-population continuous-clock
moments `f4682eb3`. It does not sample new paths or choose new extrema.

Let the two continuous uniform-order-statistic birth times be tau1<=tau2,
`C=(tau1+tau2)/2`, and `W=tau2-tau1`. All unperturbed means are full-population
means in one physical orientation. H is the derivative under the common label
policy and one of the fixed plus/minus marks.

## An exact three-part dipole

The rank-one interval gives

```
I0 = integral H(E(p)) dp       = -H(W)
I1 = integral p H(E(p)) dp     = -H(CW)
dipole_at_pref = I1-pref I0
               = -muW H(C) -(muC-pref) H(W) -H(Cov(C,W)).
```

Thus the centered first thermal moment separates into a center-displacement
term, a mean-lifetime term and a connected spread term. Their errors must be
propagated jointly; subtracting separately reported errors is invalid.

The moment derivatives are

```
H(Cov(C,W)) = H(CW)-muC H(W)-muW H(C)
H(Var C)   = H(C^2)-2 muC H(C)
H(Var W)   = H(W^2)-2 muW H(W).
```

Compute these products separately for first and second orientation. Only then
form S=(first+second)/2 and D=(first-second)/delta_cos4. In every delete-one-batch
replicate, pool the retained raw means before recomputing the products. This
preserves all shared source and mean uncertainty.

## Marginal spread versus joint fluctuation

The apparently connected C/W quantity obeys a further exact identity:

```
Cov(C,W) = (Var(tau2)-Var(tau1))/2.
```

Its derivative measures unequal broadening of the two birth marginals. It
does not, by itself, identify information in their copula. The distinct joint
direction is

```
Cov(tau1,tau2) = Var C-Var W/4.
```

Both are recoverable from the same five supplied full-population moments.
If v1=Var(tau1) and v2=Var(tau2) remain positive in all original-batch deletions,
retain the dimensionless derivative

```
H(rho12) = H(cov12)/sqrt(v1 v2)
           -rho12/2 [H(v1)/v1+H(v2)/v2].
```

This is a normalized response coordinate, not a fitted model or inferred field
count. Its source parameter remains the original common-label mark strength.

## Continuous moments

Given integer ranks k1<=k2 and n1=N+1, n2=N+2, the moment reader integrates over
the common uniform order-statistic clock:

```
E(C|ranks)   = (k1+k2)/(2 n1)
E(W|ranks)   = (k2-k1)/n1
E(CW|ranks)  = [k2(k2+1)-k1(k1+1)]/(2 n1 n2)
E(C^2|ranks) = [(k1+k2)^2+3 k1+k2]/(4 n1 n2)
E(W^2|ranks) = (k2-k1)(k2-k1+1)/(n1 n2).
```

These are not the squares of the conditional mean clock times. The reader
supplies them once; this coordinator consumes only its final committed batch
vectors and the saved thermal moments.

The formulas above were fixed in `4e708aa9` before the final baseline-moment
archive was received.

## The weak dipole contains a resolved cancellation

For **N425 minus->D**, the three dipole pieces are

| Component | Response +/- original-batch SE |
|---|---:|
| center displacement | +5.42584e-7 ± 7.55604e-8 |
| mean lifetime | −1.64931e-9 ± 1.88437e-9 |
| connected spread | −3.49536e-7 ± 9.69285e-8 |
| total dipole | +1.91398e-7 ± 1.17965e-7 |

The connected term is resolved at about 3.6 SE even though the net dipole is
weak. It is the negative of the first/second marginal variance-imbalance
response. The response of Var(C) is −1.05259e-6 ± 1.42045e-7, whereas Var(W)
changes by +1.15559e-7 ± 1.89617e-7 and is unresolved. This supports a change in
the distribution of clock centers, not just a uniform displacement of all
clocks. A prefix-dependent displacement can itself change these covariances;
the result does not eliminate that broader mechanism.

At N325 minus->D, the center term is +6.52108e-7 ± 1.66388e-7 and the connected
term −1.59355e-7 ± 1.74872e-7; the second term is not resolved there. The mean
lifetime contribution again is small. For plus->S the net dipoles are
−2.04556e-7 ± 8.26059e-8 and −2.00956e-7 ± 7.01655e-8, led by center terms
−2.53433e-7 ± 4.04909e-8 and −2.60162e-7 ± 4.15014e-8; their connected terms
remain weak.

## Joint fluctuation survives removal of timing noise

The continuous total covariance contains a conditional order-statistic term.
Writing mu1=E(tau1) and m12=E(tau1 tau2),

```
Cov_total(tau1,tau2)
  = Cov(K1,K2)/(N+1)^2 + (mu1-m12)/(N+1).
H_orderstat = [H(mu1)-H(m12)]/(N+1).
H_intrinsic_rank_cov = H_total_cov-H_orderstat.
```

The second term is the expected covariance of the shared uniform clock given
the discrete ranks. Its response is small and opposite in sign to the main
effect. The intrinsic rank-covariance response, still in continuous clock units,
is

| Source/observer | N325 | N425 |
|---|---:|---:|
| plus->S intrinsic rank covariance | +4.59145e-7 ± 7.51861e-8 | +5.34794e-7 ± 6.33625e-8 |
| minus->D intrinsic rank covariance | −1.19492e-6 ± 3.25239e-7 | −1.08880e-6 ± 1.61897e-7 |
| plus->S continuous rho12 | +6.71241e-5 ± 2.42816e-5 | +1.08497e-4 ± 2.04602e-5 |
| minus->D continuous rho12 | −2.22375e-4 ± 7.95354e-5 | −2.38442e-4 ± 5.62496e-5 |

For N425 minus->D the total covariance derivative is
−1.08148e-6 ± 1.60878e-7 and the order-statistic derivative only
+7.32127e-9 ± 3.03203e-9. Thus the joint fluctuation result is not produced by
the shared timing term. All correlation denominators are positive in both
orientations and all twenty original-batch deletions.

These signs refer to the named source/observer combinations. In particular a
negative minus->D response is a normalized geometry difference of responses,
not a statement that the unconditional covariance itself is negative. Pearson
correlation is a joint second-order coordinate, not a full copula identifier.

## Normalized plateau and one covariance ledger

The independently implemented, same-source plateau readout `1e8549b5` is
appended as supplied LOO columns, without recalculating its normalized moments.
It reports plus->S plateau variance changes +3.17944e-7 ± 1.28394e-7 and
+3.62163e-7 ± 1.04955e-7. At N425 minus->D, its variance change is
−4.41008e-7 ± 1.69308e-7 and its centroid-minus-unweighted-center response is
+8.08586e-6 ± 2.25570e-6. These describe the same shape response with positive
baseline R1 mass normalization; they do not divide by the weak derivative of
the integral and do not form independent evidence.

The result keeps the original twenty batches, full covariance with `b582015e`,
all new raw moment/thermal columns, orientation-derived LOO values, and supplied
plateau LOO values. No inverse, fit, new sampling, DP or test suite was run.
Reproduce with `scripts/p334_euler_dipole_connected_clock.py` in the managed
research Python environment. Results and factors are in
`results/p334-euler-dipole-connected-clock/`.

Scientific card: a weak E dipole can hide a resolved competition between center
movement and marginal spread imbalance. Separately, both sizes show a signed
intrinsic rank-covariance response under the same instantaneous-rank/Euler-
preserving controls. Dependency remains the original `e32a8593` fork block.
This is a finite-source response decomposition, not a new field count or an
independent replication of the preceding curve result.
