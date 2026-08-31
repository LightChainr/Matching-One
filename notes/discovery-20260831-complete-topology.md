# The complete topology readout cancels the large checkpoint-membership step

This round completes F1/F2/A_top on the original40k paired paths. The main
new numerical observation is sharp: the four R1-membership states that
explained99.85% of the restricted integrated readout explain only about
0.02%/0.005% of complete A_top's integrated orientation-contrast variance.
The full observable and the arbitrary origin allocated among checkpoint
strata make the difference; this is not a99.85% variance-reduction algorithm.

## 1. Completed source, no new random population

[9c495ab1](https://github.com/LightChainr/Matching-One/commit/9c495ab13e65f2bc93dc0849ee3b73f88724c4b1)
adds K1,K2 and all R0/R1/R2 checkpoint ranks for every original counter.
The unchanged original C++ engine took0.709s for N325 and0.963s for N425.
No geometry-pilot statistics, reliability DP, or new random samples were
computed. The20 original batches per size remain the inference units.

The complete canonical readout is `A=F1+F2-1`; its integral is exactly
`1-(K1+K2)/(N+1)`. The conditional replacement uses a common **ordered**
prefix. Both orientations must already have rank at least one, and every
required R1 clock must belong to the old whole-pair exact acceptance.
Otherwise the entire global paired vector stays original. This is the
[immediately executable global policy](https://github.com/LightChainr/Matching-One/blob/af87c4ef7ba2848e2888c78d504c48361d85092b/notes/p334-global-two-birth-loading-policy.md),
not a new R0 solver. It preserves the full covariance identity.

## 2. The completed global result

The complete source/covariance result is
[3edc785a](https://github.com/LightChainr/Matching-One/blob/3edc785a0312e4dce688bc6966593780907abc51/notes/p334-full-global-conditional-clock.md).

| Complete safe-hybrid H4-normalized contrast | N325 | N425 |
|---|---:|---:|
| A at p_ref, mean +/- original-batch SE | -.00011974 +/- .00534348 | .01010732 +/- .00694886 |
| Integrated A, mean +/- original-batch SE | .000108672 +/- .000751799 | .000943034 +/- .000819725 |
| Estimated canonical residual noise removed | 2.317% | 2.323% |
| Estimated integrated residual noise removed | 4.632% | 4.611% |

The global means are not sharply determined by this20k-per-size source.
The structural noise comparison is nevertheless informative. The same four
binary R1 flags explain99.84895%/99.86380% of the old gated integral's
individual variance, but only0.0222914%/0.0054541% of the complete hybrid
integral. This concerns that same four-state partition, not an assertion
that all nine full-rank states lack information.

All three checkpoint layers must be added with their cross-covariances.
The [exact origin identity](https://github.com/LightChainr/Matching-One/blob/467653438bfa74dbb96f988b6da5c04f07ca0f0e/notes/p334-stratum-origin-and-global-cancellation.md)
states that subtracting the common fixed alpha=1-2p_ref before stratification
leaves every complete orientation contrast unchanged. Numerically, the sum
of the three marginal stratum variances falls from.0869747 to.0122105 at
N325, and.0622048 to.00753563 at N425. The complete variances remain exactly
.0146283 and.00906587. The change is a reallocation between stratum variances
and covariances, not a reduction of the total uncertainty.

## 3. First-birth coupling changes the marked source direction

[2dd865f0](https://github.com/LightChainr/Matching-One/blob/2dd865f0b26a4d5d43f52b300293016e6ffd19b8/notes/p334-marked-global-topology-loading.md)
attaches the original direct/collective final-source mark to **full** A.
For an accepted R1 prefix and a source s,

```
A_s(p) = F2_s(p) - pi_s*[1-F1(p)]
A_direct_integral = [H2*E[T] - K1*pi_direct]/(N+1).
```

The K1 times winning-probability term is a source-weighted first-birth
coupling, not the product of separately averaged K1 and pi. It reverses both
integrated source point contrasts at both sizes on the **same accepted
pairs**:

| Size / source | Completion F2 contrast | Complete marked A contrast |
|---|---:|---:|
| N325 direct | +.0005547302 | -.0000622091 |
| N325 collective | -.0001240242 | +.0001441617 |
| N425 direct | -.0003571092 | +.0002415090 |
| N425 collective | +.0008358653 | -.0004405047 |

These are dependent archive point estimates, not established population
source signs. All unmarked R2 and complete fallback contributions remain
in a signed remainder; they are not zero or an R0-only channel. That remainder
opposes the accepted-R1 canonical sum at N325 and carries most of the N425
canonical point contrast. It is part of the full joint covariance.

## 4. Where lifetime information re-enters

With C=(K1+K2)/2 and W=K2-K1, the full integral is `1-2C/(N+1)` and cancels
W. The [next exact thermal moment](https://github.com/LightChainr/Matching-One/blob/f0dbc070b826761b097171315ab64750dee90823/notes/p334-global-thermal-moment-hierarchy.md)
recovers lifetime through a specified, non-fitted coordinate:

```
J1 = integral p*A(p) dp
   = 1/2 - [C^2+C+W^2/4]/[(N+1)(N+2)].
```

The lifetime-squared term must be separated from the center second moment.
Every higher polynomial moment uses only even powers of W, whereas the
plateau F1-F2 uses odd powers. This is a finite two-birth exchange identity,
not a field-count or asymptotic-exponent claim.

The fixed J0/J1-center/J1-width batch vectors have already been read on all
original paths and committed at
[e64febe4](https://github.com/LightChainr/Matching-One/commit/e64febe4ff10ca9cfb2f094c1b8ee8f733177fe1).
They share the global source and are not additional independent evidence.

The [shared-batch source and thermal readout](https://github.com/LightChainr/Matching-One/blob/bc05b5e41ac163649d6f9095d16e42f65b3b722a/notes/p334-global-source-thermal-joint.md)
now supplies a more specific shape direction:

| H4-normalized first thermal moment | N325 mean +/- SE | N425 mean +/- SE |
|---|---:|---:|
| Center term | .000063761 +/- .00045516 | .000448864 +/- .000511207 |
| Lifetime-squared term | .0000066243 +/- .000017950 | **-.0000364290 +/- .0000121487** |
| Total J1 | .000070386 +/- .00045046 | .000412435 +/- .000510452 |

At N425 the negative lifetime-squared coordinate is about3.00 ordinary batch
SEs from zero, while the total J1 remains unresolved because of the larger
center uncertainty. It is a raw second moment of lifetime, not lifetime
variance alone. The same-path center subtraction exposes this coordinate;
it does not create an independent replicate or establish an asymptotic
scaling law. N325 does not resolve that component at this precision.
These J1 coordinates use the original baseline paths, not a substitution
of the squared conditional mean for the required second moment.

## 5. How much clock information is in the very next direct gate?

[f3426e11](https://github.com/LightChainr/Matching-One/blob/f3426e11d1c0e67b19b1fa91b631fa9a2300590d/notes/p334-first-step-doob-clock-innovation.md)
uses the saved complete univariate clock laws. For d remaining sites,
h original direct gates, terminal readout mean m and first-step absorption
value a, the binary direct-versus-safe Doob variance is exactly
`B=h/(d-h)*(a-m)^2`. It is also a lower bound for the information in the
full next label.

On the solved real prefixes, this binary event explains6.86%-8.36% of
integrated conditional clock variance and19.33%-22.44% of the canonical
variance. The remainder requires safe-label identity and/or later suffix
information; it need not all persist after a fully observed next label.
Only already available clocks were used, including50 complete marginal
clocks inside old paired fallbacks. Coverage and missing variance mass are
kept explicit; these marginal bounds are not added into an H4 bound.

## Lifecycle and next scientific step

All numerical branches above depend on the original e81dd59 paired source.
Full births9c495ab1 and conditional clocks0d1e586d are complementary readouts
of it, not replications. The global coordinator retains shared original-batch
covariance for full observables, stratum sums, source terms and thermal
moments. No large joint covariance is inverted to manufacture an omnibus.
The final source/thermal join has56 coordinates in the same20-batch block.
Marked integral source SEs are .00069628/.00050942 (N325 direct/collective)
and .00064990/.00049139 (N425): the observed sign reversals relative to the
completion terms are not resolved population source signs.

The next route is the complete center/lifetime-shape plane and the marked
first-birth/winner coupling, rather than extrapolating the variance allocation
of a restricted R1 contribution. No new remote job, GPU rental, task closure,
PR merge or history rewrite was made in this delivery.
