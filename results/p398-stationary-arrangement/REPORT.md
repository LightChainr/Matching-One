# Dual block sizes explain most statics; arrangement still controls transitions

## Technical summary

The completed width8 stationary-response hierarchy gives a useful positive
compression and a concrete dynamical obstruction. The complete **primal and
Kreweras-dual block-size profiles explain 99.7344% of stationary-score
variance**, using60 profile classes rather than1430 labelled states. Their
remaining arrangement component is only0.2656% of score variance, but supplies
**10.17% of the final +− integrated response** because stationary and generator
contributions partly cancel.

The same60 profiles do not determine the exact transition law: **32 classes
contain different outgoing rate vectors**. An explicit pair with identical
primal/dual block sizes has **two versus four detach events** into the same
target class. This directs attention to which primal and dual blocks meet,
not simply their sizes. It does not invalidate the useful static compression
or identify a continuum field.

This is one .22415-second reuse of saved finite-state data, plus integer event
counting. No new stationary solve, rate point, width, Fréchet calculation,
Monte Carlo, server job or research test suite was used.

## A nested hierarchy separates size and arrangement information

Let `g=pi_prime/pi0` be the stationary score for the already defined generator
`F_eta=(1+eta)J+(1−eta)D`. The reference measure is the saved positive `pi0`.
For a state descriptor Z, its best squared-error score is
`g_Z=E_pi0[g|Z]`. This uses the known derivative: it diagnoses information
sufficiency and is **not an independently predicted physical measure**.

The fixed hierarchy is:

```text
block count B
  → full primal size profile lambda
  → ordered pair (lambda, lambda_dual)
  → full circular arrangement modulo rotation and reflection
```

Profiles are sorted multisets, not a list with an arbitrary block label order.
The dual profile belongs to the explicitly implemented Kreweras complement.
The total score variance is12.87669987 under `pi0`.

- **Block count:** eight classes explain89.0127%; residual10.9873%. This is the
  best arbitrary block-count function, distinct from the preceding fixed
  `−2(B−4.5)` candidate whose relative RMS error was34.27%.
- **Primal sizes:**22 classes add6.50595 percentage points, leaving4.48138%.
  Complete sizes therefore improve on count but do not exhaust the score.
- **Primal and dual sizes:**60 classes add another4.21578 percentage points,
  leaving0.265594%. This is a strong finite static compression.
- **Circular arrangement:**130 dihedral classes recover the remaining
  0.265594%; within-orbit residual is at numerical roundoff,1.24e−30 of total
  variance. Rotations/reflections are not additional stationary mechanisms.

For example, these two configurations share both profiles
`lambda=(4,2,1,1)` and `lambda_dual=(3,2,1,1,1)`:

```text
(0,0,0,1,2,1,3,0): pi0=.000184512456, score=−.058456132
(0,0,1,1,1,2,3,1): pi0=.000256769082, score=+1.949185359
```

They have opposite score signs despite identical size information on both
sides. The statewise gap is2.00764149. It is a deterministic witness selected
from this finite enumeration, not a sampling significance statistic.

## Pairing primal and dual sizes restores the source symmetry

The true rate score is Kreweras-odd. Block count is a symmetry-compatible
descriptor because complement swaps B with9−B. A primal-size profile alone
is not closed under complement; its conditional projection has maximum
oddness discrepancy4.32124 and creates nonzero diagonal normalized-response
derivatives that the full source does not have.

Conditioning on the ordered primal/dual pair fixes that problem: its maximum
oddness discrepancy is3.82e−14. This makes the paired profile a more natural
coarse source coordinate than primal sizes alone, even before considering its
better approximation. Correct symmetry still does not ensure exact closure.

## A small arrangement variance has a larger integrated effect

Keep the archived two rays, normalization and sixteen-attempt clock. For each
projected score, replace only the stationary metric derivative; retain the
same true saved generator contribution. The −+ / +− entries below refer to
archived matrix indices `(0,1)` / `(1,0)`.

The paired-size reconstruction gives integrated derivatives
`−.05351798 / −.04512655`, compared with the actual
`−.05300341 / −.05023378`. The remaining arrangement contributions are
`+.0005145631 / −.0051072292`: respectively .97% and10.17% of the magnitudes
of the final responses. Both reconstructions now have the correct signs;
the preceding block-count-only reconstruction did not.

The +− response has this additive attribution, rounded here for readability:

```text
generator contribution          +.06552032
block-count score                −.05256505
additional primal-size score     −.03502478
additional dual-size score       −.02305704
remaining arrangement score      −.00510723
total                            −.05023378
```

These are contributions to one response, not independent pieces of evidence.
The score increments are orthogonal in `L2(pi0)`; their images in the selected
response need not be orthogonal or proportional to their score variances.
Cancellation explains why a small omitted score component can matter to the
readout. Additive reconstruction agrees to3.86e−16 in the matrix norm used
by the script. Static cross derivative also improves: paired sizes give
.78684149 versus the full .78211719.

## Two versus four physical detach channels

An exact Markov compression requires every microstate inside a profile class
to have the same total transition rate into each other profile class. Integer
event counts show failure in32 of the60 classes, both for `J+D` and `J−D`.
For the following pair, the source profiles are identical:

```text
x=(0,1,0,2,3,4,2,5)
y=(0,1,2,0,3,4,5,3)
lambda=(2,2,1,1,1,1), lambda_dual=(3,3,2)
```

Into the target `((2,1,1,1,1,1,1),(5,3))`, x has two detach events and y
has four; neither has a join event into that target. Consequently the exact
rates are `2(1−eta)` and `4(1−eta)` in the inherited attempt clock.
This is an off-diagonal target class, so no diagonal convention is involved.

The difference has a simple incidence interpretation. Connect each site to
its primal block and its dual block. In x, the two primal size-two blocks meet
dual sizes `(2,3)` and `(3,3)` respectively. In y, both meet `(2,3)`.
Splitting either endpoint of a primal pair meeting dual sizes2 and3 merges
those dual blocks into size5. Thus x has one eligible pair and y has two,
giving the observed two-versus-four channel count. Equal degree/size lists
do not say which degrees meet.

The excluded object is the **exact60-profile Markov quotient for arbitrary
initial microstates**. This does not exclude an approximate profile model,
special stationary coarse behavior, a useful memory kernel, or a compact
model with explicit incidence information.

## Definitions, provenance and interpretation boundary

Inputs are the [completed rate-response package](https://github.com/LightChainr/Matching-One/blob/fb01c44aa45e4f8d37d52144e2ad7c4adfe6ce40/experiments/p398-linear-response-20260831/README.md)
at `fb01c44a` and the preceding
[block-count score archive](../p398-block-count-measure/REPORT.md) at
`7da1eeb0`. The unchanged constructor supplies state ordering, complement,
basis and exact join/detach actions. Saved1430-state `pi0`/score and the
186-dimensional character-i matrices supply the floating-point response.

For each projected score h, set `Bprime=Q* diag(pi0 h) Q`. With saved rays Z,
mass H and normalized integrated response R, the stationary derivative is

```text
C0prime(h) = conjugate[Z* Bprime Z]
Rprime_pi(h) = C0^-1 {conjugate[Z* Bprime H^-1 Z] − C0prime(h) R}
```

The true saved generator derivative is then added, unchanged. This separates
measure information from already computed propagation; it does not construct
a different autonomous stochastic model. The variance decomposition and
stationary scores are deterministic float64, while the event-count witness
is integer arithmetic. No rational certification of stationary probabilities,
Monte Carlo uncertainty, finite-strength law, site-Matching operator map or
continuum-field count is claimed.

The [contract](../../analysis/p398_stationary_arrangement_contract.json) and
[script](../../scripts/analyze_p398_stationary_arrangement.py) were committed
before calculation at `c643249d01abb402d377f6b319a66a1433aeca65`.
[The result JSON](latest.json) records input/code hashes, all state descriptors,
all conditional scores, exact witnesses and response components. To reproduce,
use a separate worktree at that execution commit, which predates the saved
output, ensure the pinned #509 source objects are fetched, and run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python scripts/analyze_p398_stationary_arrangement.py
```

## Next mechanism: which size classes meet?

The specific next candidate is **size-conditioned primal–dual incidence**,
motivated by the exact two-versus-four witness. A count of size-two primal
blocks joining dual sizes2 and3 predicts this particular transition rate;
that local success is not yet a complete stationary-score law. A useful next
calculation would ask whether such defined incidences predict the omitted
arrangement response, preferably with a coefficient or transition identity
derived from the generator instead of another unrestricted score fit.

This finite-model line remains subordinate to the original norm4 mechanism
question. Translating it requires an explicit microscopic source/readout map.
P334's independently delivered contact-loading results offer a related question
about which components meet, but they concern different objects and samples;
their resemblance is a hypothesis bridge, not shared empirical confirmation.
