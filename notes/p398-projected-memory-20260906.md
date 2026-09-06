# #588 Phase A — the held-out failure is not one defect, it is two

`scripts/p398_projected_memory.py`, `results/p398-projected-memory/latest.json`,
`tests/test_p398_projected_memory.py` (19 tests).
Frozen input: #580's rank-6 observable-Krylov span at `eta = 0`, quoted through the
same code path rather than rebuilt.

---

## The question

#580 left a fork its own language could not close. A frozen rank-6 span transports a
declared intervention for `+0.006` excess relative error from width 4 to width 8, and
fails badly on five readouts it was not built from — already at `eta = 0`. Two
completely different things produce that pattern:

```text
(i)  the dictionary is missing an instantaneous state coordinate;
(ii) an exact Markov process, projected onto an observable-relative subspace,
     generates memory.
```

They call for opposite repairs, and everything downstream — #584's index, #249's
Hankel rank, any eventual square-site reading — inherits whichever is true.

P398 settles it because the microscopic generator is known exactly, so the memory
kernel is computed rather than estimated. With `P = Phi Phi^T`, `Q = I - P`, and `G`
acting on functions:

```text
A = Phi^T G Phi     B y = Phi^T G y  (y in range Q)
C a = Q G Phi a     D   = Q G Q

d x_R/dt = A x_R(t) + int_0^t K(t-s) x_R(s) ds + B exp(tD) x_U(0)
K(tau)   = B exp(tau D) C
```

`K` is exactly the discarded dynamics. The third term is exactly the part of the
readout the span does not see at `t = 0`. **They are different defects and this run
separates them.**

## What had to be true before any number counted

`D` is 1424 x 1424 at width 8 and is never formed: the columns `c_j = Q G phi_j` are
propagated as full-space vectors constrained to `range(Q)`, `exp(tau D)` is taken by
uniformization on `I + D/Lambda`, and `K` is read off by two inner products. That is
fast, and it is also exactly the kind of shortcut that returns something smooth,
plausible and wrong. So:

```text
full-rank projection, K must be identically zero      0.000e+00   (exactly)
generator-closed Krylov span (dim 10 of 14), K = 0    6.098e-30   vs ||A|| = 19.05
synthetic planted (A,B,C,D), no vanishing row sums    3.601e-16
explicitly built D versus the matrix-free path        3.177e-14
propagation leakage back into the span                5.797e-16
resolvent Schur identity, width 4                     7.427e-16
resolvent Schur identity, width 5                     6.276e-16
```

The resolvent check is the one that matters most: `Phi^T (sI-G)^{-1} Phi` must equal
`(sI - A - B(sI-D)^{-1}C)^{-1}`. Transposing `B` and `C` leaves every norm in this
file unchanged and still yields a decaying kernel; it does not survive that identity,
and the identity contains no time integration, so a failure could not be blamed on a
quadrature rule.

## Result 1 — the memory is short, weak, and low-order, and stays that way

Declared intervention, `eta = 0`, #580's frozen rank-6 span:

```text
width  states   ||A||   ||K(0)||   int||K||   int||K||/||A||   decay   tail>2   poles
  4      14     12.99     5.47       0.876        0.067        0.063   0.0000     2
  5      42     12.53     5.81       1.016        0.081        0.089   0.0000     3
  6     132     12.35     6.61       1.207        0.098        0.105   0.0001     4
  7     429     12.40     7.67       1.456        0.117        0.119   0.0003     4
  8    1430     12.56     8.45       1.681        0.134        0.136   0.0006     4
```

`poles` is the block-Hankel order carrying 99.9% of the kernel energy; at the 99%
level it is 3 at every width from 5 up.

Two readings, and the second is the important one.

The memory **decays before the first declared lag**. Its centroid is at `t = 0.14` at
width 8; the declared lag grid starts at 0.25 and runs to 4.0, and the mass beyond
`t = 2` is `6e-4`. So the projected process is nearly Markovian *in time* — but with a
generator that is not `A`.

And the memory **does not become more complicated as the state space grows**. From 14
states to 1430 — a factor of 102 — the effective order goes 2, 3, 4, 4, 4 and then
stops. The integrated weight roughly doubles and the decay time roughly doubles, but
the pole count saturates. This is a clean negative for the "growing memory /
predictive noncompression" branch of #588's table: on this object, projecting onto a
6-dimensional observable dictionary does **not** export unbounded complexity into
history.

## Result 2 — for the declared dictionary, the entire residual is memory

The four-closure attribution is the sharpest thing in the run. The projected equation
is an identity, so with the exact `K` and the exact inhomogeneous term the reduced
solve must reproduce `Phi^T exp(tG) f`. Turning each term on alone therefore
attributes the Markov closure error exactly, with no norm convention in it. All four
use the same integrator, so the differences are the terms and nothing else.

Declared readouts, `eta = 0`, step `h = 1/16`:

```text
width    markov     +memory    +forcing      +both     memory share
  4    2.475e-03   4.612e-04   2.475e-03   4.612e-04       81.4%
  5    5.398e-03   4.613e-04   5.398e-03   4.613e-04       91.5%
  6    8.257e-03   4.633e-04   8.257e-03   4.633e-04       94.4%
  7    1.083e-02   4.683e-04   1.083e-02   4.683e-04       95.7%
  8    1.309e-02   4.744e-04   1.309e-02   4.744e-04       96.4%
```

`+forcing` equals `markov` to the digit because every declared readout seeds the
Krylov span, so `Q f = 0` and the inhomogeneous term is identically zero. That is not
a coincidence to be reported as a finding; it is a structural property of the frozen
configuration, and it is why the two blocks have to be scored separately.

At widths 4-6 the same rows at `h = 1/32` give memory shares of 95.1 / 97.8 / 98.6%
and a `+both` residual four times smaller (1.16e-04 against 4.61e-04) — second-order
convergence, which is what says the leftover is quadrature and not a wrong kernel.

**So #580's frozen model, on its own dictionary, has a residual that is 96% projected
memory, of effective order 3-4, decaying in 0.14 time units.** A memory closure with
three or four poles would essentially close the declared dictionary exactly. That is
the strongest form of #588's "few memory poles" branch, and it is Phase C's job to
find out whether those poles can be Markovianized by a few common auxiliary
coordinates.

## Result 3 — for the held-out readouts, neither repair is the object

Held-out readouts, `eta = 0`, `h = 1/16`:

```text
width    markov     +memory    +forcing      +both    mem%   frc%   resid%
  4    1.372e-02   7.786e-03   8.911e-03   7.825e-04  43.3   35.1    1.4
  5    2.489e-02   1.283e-02   1.595e-02   7.723e-04  48.5   35.9    0.8
  6    3.994e-02   2.353e-02   2.270e-02   8.671e-04  41.1   43.2    0.5
  7    5.169e-02   3.055e-02   2.923e-02   1.063e-03  40.9   43.4    2.1
  8    6.182e-02   3.768e-02   3.528e-02   1.227e-03  39.0   42.9    2.0
```

Memory alone removes about 40%. The unrepresented readout alone removes about 43%.
Together they remove 98%. Neither dominates, at any width, and the split is stable.

**The held-out failure is not "missing state" and it is not "hidden history". It is
both, in roughly equal measure, and they interact.** Adding instantaneous coordinates
would leave 40% of the error standing; adding a memory kernel would leave 43%.

This is a branch #588's decision table does not have, and the script says so in the
verdict string rather than rounding to the nearest listed one:
`MEMORY_AND_UNREPRESENTED_READOUT_ARE_JOINTLY_REQUIRED`. The same thing happened in
#580, whose table lacked `STATE_TRANSPORTS_ON_ITS_OWN_READOUTS_ONLY`. Twice now the
honest reading has been a branch the design did not anticipate; that is worth noticing
about the design, not just about the object.

### Why the norm ratio would have said the opposite

`||int K x_R|| / ||A x_R||` is 0.029 for the declared readouts and 0.062 for the
held-out ones at width 8 — both tiny. Read alone, those numbers say memory is
negligible, and they are the numbers the analysis plan declared as the primary
discriminator. They are misleading because `A` carries the fast modes, so the ratio
measures memory against the largest term in the equation rather than against the
*error*. The attribution replaces it, and the plan's thresholds are reported alongside
rather than quietly dropped.

Provenance, GOVERNANCE §2C: the norm-ratio thresholds were fixed in
`notes/analysis-plan-20260906.md` before any kernel existed. The four-closure
attribution and its 0.5 shares were added after the width 4 and 5 controls were read
and before widths 6, 7 and 8 were computed. 0.5 is the natural half, not a tuned value,
and the verdict does not sit near it: 0.39 and 0.43 are both clearly below.

## Result 4 — the memory is not pencil-sensitive; #580's positive realization is

A3 repeats everything for `single_point_join`, which tilts one boundary point and lies
outside the `span{J,D}` pencil that built the state space.

The raw numbers look like a large pencil effect and are not one. The two interventions
are not the same size: at equal `eta`, `uniform_join_minus_detach` tilts all `2w`
moves and `single_point_join` tilts one, so `||H||/||G_0||` is 0.396 against 0.099 at
width 8. Per unit of perturbation:

```text
width   in pencil   out of pencil
  4       2.465        1.325
  5       1.878        1.678
  6       2.228        1.803
  7       2.429        1.810
  8       2.577        1.799
```

A factor of 1.4, not the factor of 5.7 the raw kernel distances show. **The memory
kernel responds to an out-of-pencil intervention at essentially the same rate as to
the declared one.** The baseline memory subspace still carries the intervened kernel
in both cases (residual 0.068 in pencil, 0.011 out of pencil, at width 8), and the
pole count is unchanged at 4.

Set that against #580's positive-realization result on the identical pair: the
coarsest exact strong lumping valid for the whole affine family is 750 blocks in
pencil at width 8, and **collapses to the identity partition — all 1430 blocks — out
of pencil, at every width.** That is all-or-nothing and cannot be explained by
perturbation size.

So the three descriptions of the same object respond differently to the same probe:

```text
r_transport   (signed low-rank, frozen)   survives the out-of-pencil probe (#580)
r_positive    (exact lumping)             destroyed by it (#580)
memory        (order, subspace, decay)    barely notices it (#588)
```

Reporting one number called "state dimension" for this object would be reporting a
choice of which of these three to look at.

## What this does and does not settle

Settles:

- The declared dictionary's residual is memory, it is low-order, and it is short.
  A small state plus three or four poles is a genuine candidate description.
- Memory complexity does not grow with width on this object. The predictive
  noncompression branch does not fire here.
- The held-out failure is jointly memory and unrepresented readout, ~40/43, stably
  across a factor of 100 in state count.
- The memory description is not an artifact of the operator pencil.

Does not settle, and is Phase B/C:

- Whether the 3-4 poles can be Markovianized by a few *common* auxiliary coordinates,
  or whether they are readout-specific. That is the difference between an operational
  finite fiber and a non-Markov reduced description, and it is Phase C.
- The nested-dictionary growth profile `r_linear / E_markov / M_memory / r_memory /
  r_transport` for D0/D1/D2, which is Phase B. D0/D1/D2 membership is closed as of
  #580; no observable may be added after a kernel has been looked at.
- The owner's balanced-reduction adversary: finite-horizon Gramians on the zero-sum
  contrast subspace, stationary mode removed explicitly, balancing transform frozen at
  `eta = 0`. Phase C, and it is the strongest deterministic Markov competitor.

### For #584

A structured residual does not imply an omitted instantaneous coordinate — but on this
object it does not imply unbounded memory either. If #584 finds a plausible fiber
label, the alternative it has to beat is now quantified: a 3-4 pole, `tau ≈ 0.14`
memory kernel. That is a small, concrete competitor, not a vague warning.

## Claim boundary

Mori-Zwanzig memory is a property of a chosen projection and dictionary. It is not a
physical field. A low-order kernel is not evidence for an LCFT Jordan block. An
auxiliary Markovianization coordinate would be an input/output realization coordinate
until independent representation/modulus/topology fingerprints identify it. P398 is a
calibration model: nothing here transports to square-site Matching One without a
declared intervention with common observables before and after, and a supplied map
between microscopic state spaces.
