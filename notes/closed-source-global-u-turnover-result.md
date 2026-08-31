# The closed attractive source amplifies U weakly, then suppresses it

**Subsequent result:** the
[exact strong-coupling series](closed-source-angular-strong-coupling-law.md)
now fixes the eventual sign: U approaches0 from below. Since the last
frozen point below is positive, a later zero crossing and negative valley
are necessary. No extra point was added to this completed four-point grid.

**New finite-model conclusion:** the already selected source
`S*=C+F+Bv` does not monotonically amplify the original global U.
The exact N25 pair has `U_t(0)=+0.126165363414`, but
`U_t(log 2)=-1.370778221631`. There is therefore at least one local
maximum in `0<t<log 2`. This decision used four couplings fixed before
the new histogram or scores; no source or amplitude was fitted.

The [full score](../results/p337-closed-source-finite-coupling/score/REPORT.md)
and [rational enclosures](../results/p337-closed-source-finite-coupling/score/score.json)
are complete. The zero-coupling line is imported, not recalculated.

| exp(t) | critical p | U | U_t |
|---:|---:|---:|---:|
| 1, previously published | .592665539328 | .880466156963 | +.126165363414 |
| 2 | .710261724749 | .242368571338 | -1.370778221631 |
| 4 | .809483564856 | .000283785592447 | -.00340200365326 |
| 8 | .890410580736 | 4.58979443475e-8 | -5.97094985292e-7 |
| 16 | .941391938090 | 3.97264799852e-12 | -5.70475721983e-11 |

All four new U_t/A enclosures lie strictly below zero. This bounds one
local peak's location, not its height, uniqueness or the location of
the global maximum. We stop at the frozen grid.

## Why attraction can extinguish the oriented signal

The exact action can be written

```text
J = beta1 + beta_null = 2 beta1-r,
S* = J-3K+2N+1,
g = 2K-J = Bmix-2CB+r,
weight(A) proportional to h^K m^(-g),
m=exp(t), h=p/((1-p)m).
```

Both cycle dimensions are increasing and supermodular. Their sum gives
the [positive-association theorem](closed-source-attractive-cycle-gas.md).
Together with a positive-probability pivotal configuration, this proves
the [unique simple matching root for all finite t>=0](closed-source-critical-root-order.md).
The loss of U is not a zero-slope singularity or a switch between roots.
Positive association orders q, but not E=q^2 or its signed projection.

The [strong-coupling proof](closed-source-two-state-turnover.md) uses a
cut bound, not a fit. Empty and full states have g=0; every other state
has g>=2, and rank-one states have g>=3. At fixed finite N, the critical
law approaches an equal mixture of empty and full states:

```text
h0 -> 1,  logit(p0)=t+O(exp(-2t)),
mean occupation -> 1/2,
P(rank=1)=O(exp(-3t)),
mean(E)=1-P(rank=1) -> 1,   Q_h -> N/2,
U(t)=O(exp(-3t)) -> 0.
```

The middle rank sector carrying E variation disappears while Q retains
a nonzero thermal slope. The Bernoulli coordinate p0 is not the actual
occupation density of the interacting measure. The exponent3 is a decay
bound, not an identified leading power or continuum exponent.

Thus the same fixed source has a positive weak response but removes the
rank-one population needed by the global readout at strong coupling.
The limiting argument and old positive signs force a finite positive
global maximum; the new finite grid also locates a local maximum before
log2. The bound was proposed in the freeze; the full proof was completed
in parallel without consulting the new finite-coupling values.

## A positive microscopic model and its endpoint consequences

The [two-current representation](closed-source-two-current-representation.md)
constructs this action for every integer m>=2. One circulation is arbitrary;
a second has zero total ambient winding. Their cardinalities are
m^beta1 and m^beta_null. The latter is a face-height boundary, so the model
has local nonnegative constraints. Both currents use the same occupied
graph; unconditional independence is not asserted. The site activity in
this representation is y/m^3, distinct from the two-state h=y/m.

The exact saturated identity remains
`U_parent,end(t)=2^(13/8) U_child(t)`. The unique-root proof removes its
earlier provisional simple-root qualifier for finite positive t. Nested
saturated endpoints inherit this turnover, multiplied by the fixed area
factor. These are transported predictions, not new parent measurements.

The preceding [single-defect result](checkerboard-single-defect-global-u-result.md)
excluded source-independent geometric gain away from that endpoint.
The targets are distinct: endpoint closure holds; its scalar-gain interior
extension fails; monotone source amplification also fails on the homogeneous
child. No added source is introduced to rescue either failed candidate.

## Scientific card

- **Mechanism changed:** monotone global-U amplification by this closed
  source is excluded on the fixed pair. Rank-one depletion yields U->0.
- **Observer/sector/source/geometry:** original separately normalized,
  pooled-root/slope U; ordinary q/E, normalized cos4; exp(t S*);
  N25 (5,0)/(4,3), with different Smith classes.
- **Lifecycle:** prediction/grid freeze b70dc4bd; scorer2b3b375d;
  producer693d0d4a; counts e49be207; one score a70eeff0; branch delivery.
  Companion proofs are theory, not separate random evidence.
- **Dependency:** same N25 exhaustive populations as the earlier first-
  moment packet, now with the finite-source sufficient histogram.
  Zero new random samples; no independent production claim.
- **Cost:** two 2^25 passes, 3.066seconds including compilation;
  all four scores1.508seconds. Local CPU only; no cloud or test suite.
- **Boundary:** fixed finite volume; no thermodynamic transition,
  continuum identity, fitted asymptotic power, unique peak or larger-N
  conclusion. The stopped F4 and lag-one blocks remain stopped.
- **Size follow-up:** the
  [winding-barrier theorem](closed-source-winding-barrier.md) and
  [size/sign/transmission synthesis](closed-source-size-sign-and-transmission.md)
  now supply the next size predictions and the separate fixed defect
  decision. Extra points around this resolved peak remain out of scope.
