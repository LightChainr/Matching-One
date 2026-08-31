# Shared next positions align the two birth responses and partly cancel in E

The new conditional trajectories expose a clear finite-source mechanism:
first-birth and completion responses to the **same next label** have positive
covariance. Their common mode reinforces the odd topology observer A and is
partly removed by the even observer E. This is a measured continuous-response
coupling on the paired geometries, not merely a runtime or precision metric.

## New acquisition, unchanged original prefix population

[Producer a3249a59](https://github.com/LightChainr/Matching-One/commit/a3249a598e7af6db940fd602df31533ab5c13d38)
retains all original20,000 ordered prefixes at each of N325/N425. For each
prefix it draws eight quartets: independent uniform next labels U,V (allowing
equality), and two independent remaining suffixes under each label. Each tail
is shared by the two original orientations. This produces32 new paired tails
per prefix, **1,280,000 total**, without another network DP or geometry change.

[Complete data e32a8593](https://github.com/LightChainr/Matching-One/tree/e32a85939279b8574278024d647b56d2d1485247/results/p334-nested-next-label-forks)
contains40 compressed original-batch files and both metadata records, about14MB.
N325 took16.77s and N42523.99s locally. These are new auxiliary conditional
random trajectories on old prefixes, not1.28 million independent new prefixes
or an independent replication of the old population signal.

## 1. A positive cross-birth response with an observable cancellation size

For a fixed original prefix Z and a common next label U, let m_i(U) be the
conditional paired H4-normalized response of F_i. Define

```
B_ij = E_Z Cov_U(m_i(U),m_j(U) | Z),
Gamma_next = B_12 = [B_AA-B_EE]/4,
eta_next = 2*B_12/(B_11+B_22).
```

All these B terms concern the **next-label** component of conditional
variation, not the entire population variance. The exact A/E identities give

```
B_AA = (B_11+B_22)*(1+eta_next),
B_EE = (B_11+B_22)*(1-eta_next).
```

The new data yield:

| Endpoint | N325 Gamma_next +/- SE | N425 Gamma_next +/- SE |
|---|---:|---:|
| canonical p_ref | .001381329 +/- .000044623 | .000997369 +/- .000040007 |
| p integral | .0000226478 +/- .00000150920 | .0000119935 +/- .000000823242 |

The corresponding cancellation proportions, relative to the sum of the two
separate birth-response variances, are:

| Endpoint | N325 eta_next +/- SE | N425 eta_next +/- SE |
|---|---:|---:|
| canonical p_ref | .231366 +/- .00706190 | .247563 +/- .00870600 |
| p integral | .305290 +/- .0184508 | .322626 +/- .0210751 |

Thus E removes roughly23–25% of this canonical next-site variance scale and
31–32% of the integrated scale; A receives the corresponding enhancement.
This is not a claim of a23–32% reduction in all E uncertainty or runtime.

Every nonzero contribution is in the five joint checkpoint cells containing
R0. When both orientations have rank at least one, their first births are
already prefix-measurable. The paired first-birth innovation is then exactly
zero and every A/E tail difference lies on `(1,1)`. R0 activates the missing
first-birth direction `(1,-1)` and its coupling to completion.

### Positivity is not implied by birth ordering alone

The inequality K1<=K2 and monotonicity of nested birth events do not by
themselves force positive next-label covariance. A four-label abstract example
illustrates this without an enumeration: birth1 occurs when a or b has been
inserted; birth2 when b,c,d have all been inserted. The events are monotone and
nested. Conditional on the first label a,b,c,d, their mean clock pairs are

```
(1,4), (1,11/3), (7/3,11/3), (7/3,11/3).
```

Their covariance across the four uniform first labels is **-1/18**. Integrated
F1/F2 covariance is consequently -1/450. This is only an abstract nested-event
counterexample, not a claimed torus realization or a negative H4 example.
It identifies exactly which weak assumptions do not explain the measured
positive covariance. A stronger torus-specific explanation remains a
scientific question.

## 2. Most continuation variation still lies after the next position

On the original01+10 prefixes, the next-label fraction of first-birth suffix
variance is14.54% +/- .51 percentage points at N325 and11.39% +/- .57 at N425.
For complete integrated E it is only5.128% +/- .489 and4.064% +/- .479.
These cells carry about45.9%/45.8% of all first-birth suffix variance and
33.1%/33.3% of complete E's integral suffix variance.

The whole-population next-label shares are approximately14.23%/11.60% for K1,
20.66%/17.43% for canonical E and6.98%/5.41% for integrated E. The coherent
next-site coupling is clear, but a single insertion does not exhaust the
remaining random degrees of freedom. This supports following safe-site
reorganization and later continuation, not declaring a one-gate closure.

The N425 complete E integral averaged over the32 new tails is
`-.000722835 +/- .00031093`, with the same original20-batch inference units.
This remains a negative point direction; the paired dependence on the old
baseline and exact-clock hybrid is retained, not treated as fresh prefix
evidence.

## 3. Immediate gate co-promotion is a different spatial projection

In01/10, name the two local roles R0-first-birth and R1-second-birth. For a
prefix with d vacant labels, their trigger sets G0,G1 give the exact
conditional same-position covariance

```
Cov_U(g0,g1 | Z) = |G0 intersect G1|/d - |G0|*|G1|/d^2.
```

This is excess overlap relative to the product of gate fractions **inside
each prefix**. Replacing it with a product of population gate rates would
mix in prefix heterogeneity. The new two-label half-difference estimator
reads it without enumerating all d positions.

Its full-population weighted01+10 mean is
`.00038125 +/- .00009521` at N325 and
`.000059375 +/- .00007820` at N425. Only N325 resolves positive average
co-promotion in this immediate Boolean projection. The much clearer Gamma
at both sizes uses continuous future-response functions, not just whether
the very next insertion already completes a birth.

With the fixed **unscaled low-rank minus high-rank** orientation convention,
both gate roles covary negatively with E, while their A signs oppose. But
restoring the original first-minus-second/H4 convention substantially
cancels01 against10: all four resulting gate-by-E population covariances
are below two batch SEs. Local gate co-promotion is therefore not already
an explanation of the global E mean signal.

The first-birth gate includes direct0->2 jumps. Their observed counts among
the sampled next labels are26/38,496 and26/39,312 for N32501/10, and
16/40,800 and14/38,432 for N425. Each independent next-label draw is counted
once; its two suffix rows are not extra gate observations.

## 4. An exact contact-generator description of first birth

The new [rank-zero contact theorem](https://github.com/LightChainr/Matching-One/blob/e67d9b900c535bff489c149cf0bd559acddb08c7/notes/p334-r0-singleton-contact-birth-generators.md)
identifies the winding information an insertion actually introduces. In an
ambient-rank-zero occupied component C, let p_C(u) be a root-to-contact lift
potential. For a vacant site v and its contact edge e=(v,u), set
`alpha_e=delta_e-p_C(u)`, where delta_e is the edge displacement. After choosing
one anchor contact e_C in each component, the new homology image is exactly

```
span_Z { P^{-1}(alpha_e-alpha_e_C) : e contacts the same component C }.
```

Differences are taken **within** components; subtracting two independently
gauged component addresses would invent a winding generator. Contracting old
spanning trees and the anchor star leaves precisely these fundamental cycles.
The construction uses the existing union-find displacement information, not a
new exhaustive solver.

For the present four-distinct-neighbor NN geometries, a direct R0->R2 jump
has two possible contact architectures: one component supplies at least three
contacts with two independent winding differences, or two components supply
two contacts each with independent winding lines. The note gives explicit
5x5 examples of both, as well as a three-contact example with zero ambient
birth. Thus contact count and component partition alone cannot classify the
birth: the relative lift addresses matter.

These are exact structural results and paper examples. The rare direct jumps
counted above have **not yet** been classified into these architectures. The
new bridge to the positive Gamma is to attach these inexpensive contact
coordinates to next-site responses, including sites with no immediate birth.

## 5. Finite fork averages retain a specific amount of conditional noise

For one quartet put a=X_U0-X_V0 and b=X_U1-X_V1. The recorded matrices are

```
Vhat=(aa'+bb')/4,
Bhat=(ab'+ba')/4,
What=(a-b)(a-b)'/4,
Vhat=Bhat+What.
```

They estimate total suffix covariance V, next-label covariance B, and the
remaining covariance W after observing that label. Bhat can be indefinite
on a finite sample and is not clipped. With16 independent next-label groups
and two suffixes per group, the32-tail average retains `B/16+W/32`, not V/32.
The actual variance gain relative to a single old tail is
`15B/16+31W/32`. The old-minus-fork residual instead has second moment
`17B/16+33W/32`; calling that a gain would be incorrect.

These [exact nested-fork identities](https://github.com/LightChainr/Matching-One/blob/84018e1969ce6cea9537a21f26f01b155f4f3afd/notes/p334-nested-label-fork-innovations.md)
preserve the common ordered prefix and the original covariance units.

## Reusable scientific outputs

- [Complete nine-cell Doob analysis and batch covariance](https://github.com/LightChainr/Matching-One/blob/24872eef688f2aad2c1288ccbc74a3928aacbc30/notes/p334-next-label-doob-result.md).
- [Gamma, eta and shared gate covariance](https://github.com/LightChainr/Matching-One/blob/e0494fdf1b2cefdf899a5663e6d60477257d3ab5/notes/p334-next-label-eta-gate-joint.md).
- [Immediate-gate spatial projection and direct-jump counts](https://github.com/LightChainr/Matching-One/blob/c6ee37a84eb99a68513e4a81f6aa67870e131b86/notes/p334-next-label-gate-coupling.md).
- [Exact nested-fork theory](https://github.com/LightChainr/Matching-One/blob/84018e1969ce6cea9537a21f26f01b155f4f3afd/notes/p334-nested-label-fork-innovations.md)
  and [rank-zero contact generators](https://github.com/LightChainr/Matching-One/blob/e67d9b900c535bff489c149cf0bd559acddb08c7/notes/p334-r0-singleton-contact-birth-generators.md).

## Lifecycle and next discriminant

The population source is still e81dd59f/9c495ab1, now with fresh auxiliary
suffix RNG at e32a8593. All complete nine-cell contributions use full20k
denominators. One covariance coordinator retains shared original-batch
uncertainty for mean responses, full noise matrices, signed Gamma, eta and
the26 gate coordinates. No high-dimensional covariance inverse or PSD clipping
is used. Immediate-gate alignment depends on the fixed common-label coupling;
no intrinsic spatial covariance, asymptotic exponent or new field identity
is established by these results.

The next discriminant is whether the continuous positive Gamma persists among
shared next labels that leave both checkpoint ranks unchanged, and which
rank-zero winding-contact patterns carry it. That distinguishes reorganization
before a visible birth from immediate gate effects, using the saved new forks
before any further trajectory expansion.
