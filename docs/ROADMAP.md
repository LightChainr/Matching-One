# Roadmap

This roadmap optimizes for **information gained per unit effort**. It is not a permission system. Existing-data analysis, exact work, pilots and exploratory production may proceed whenever useful.

The only hard constraints are chronology, observable-semantic compatibility for claim-bearing scores, and non-duplication of correlated evidence.

For a fast lifecycle view—what is already implemented, what is actually missing, and the cheapest next information—see `docs/NEXT-TARGETS.md`.

## Default attention now

The useful question is not “which issue number is P0?” but “which next result would most change the live model space?”

### 1. Resolve the lower spin-4 representation/selection question — #257/#262/#250/#264

The continuum preflight is already complete:

```text
V_(2,+/-2): x=17/4, |spin|=4, four legs
Q4 epsilon:  x=21/4, |spin|=4
```

and their `Q` velocities are distinct. What is missing is the physical `Q=1` Potts projector/multiplicity and the actual lattice overlap.

**Best next information:** an exact zero/nonzero global matching matrix element, projector/multiplicity result, or controlled charged four-leg insertion. This is higher value than another free radial exponent fit because it can remove or force accommodation of the lower-dimensional competitor.

### 2. Calibrate the logarithmic diagnostics on the known triangular pair — #246/#234

The exact protocol, tiny oracle and production implementation already exist. The remaining value is a small scaling calculation on a logarithmic pair whose lattice definition comes from outside the square-site interpretation.

**Best next information:** show whether the same rank/Jordan diagnostics that are being proposed for square-site data recover the known energy/log pair. A failure here should improve the diagnostic before it is used to name the unknown target sector.

### 3. Same-N norm-5 coalescence — #205

The exact A/B/C geometry and frozen H4/H8/H12 residuals are ready. At fixed N the new C node changes Smith class while the primary H4 interpolation removes radial exponent, parent amplitude and `p_c`.

```text
N325: 5 M_C - 11 M_A + 6 M_B = 0
N425: 20 M_C + 13 M_A - 33 M_B = 0
```

**Best next information:** a variance-only pilot followed by the smallest fresh common-field target that can decide whether H4 interpolation survives the arithmetic-class change.

### 4. Implement the minimal-state analysis on existing blocks — #249/#180/#253/#255

The repository already has low-rank matrix/semigroup and structured Krawtchouk assets, but the context-Hankel/minimal-realization program itself is not yet implemented on `main`.

**Best next information:** covariance-aware held-out rank/composition on archived norm-2, norm-5, N290 and compatible local/score-mode data. This is no-new-production work and can tell us whether several current descriptions are genuinely distinct mechanisms or different coordinates of one small predictive block.

### 5. Turn the modular-scalar theorem into one typed square-site score

The new exact modular channel classification says that rank/either/cross are modular scalars, while direction-specific channels are basis-dependent. This opens a cleaner scalar elliptic-spin route.

Keep this separate from the older primitive C3/Pell experiment: that block found a real primitive character but **failed** the simple ordinary-H4 zero bridge.

**Best next information:** one typed `cross` or `either` square-site observable at a new modulus with a parameter-free Q4/Jordan-versus-competitor shape prediction.

### 6. Re-rank norm-4 before buying the full production — #154

The N260/N340 variance pilot and scorer are ready. Under the current source/target plan the expected q2/Jordan separation is useful, but the production is expensive and source covariance dominates.

**Best next information before compute:** rerun the already-existing information/CPU optimizer with today’s enlarged model set. Norm-4 is ready to run; readiness alone does not make it the best next purchase.

## Parallel theory / exact lines

These do not need to wait for the default attention order.

### Matching observable -> continuum bridge — #61/#114/#233/#120

There is already substantial structure:

- exact finite matching/complement semantics;
- exact Russo/pivotal derivative relation;
- an explicit C4 self-matching microscopic tangent with exact UV involution;
- FK/Potts homology and thermal-cumulant formulas;
- transfer/defect formulations proposed as stronger bridges.

The missing step is not another statement that the sectors look even/odd. It is an RG/FK/defect/transfer-matrix map—or a precise obstruction—that predicts a field-selection rule.

### Metric-free cross-microscopic ratios — #118

The dimensionless combinations are already derived. The next information is empirical universality across a second microscopic realization at matched modulus and typed observable semantics.

### Threshold-origin line — #123/#3/#13/#17/#29/#47

Keep this scientifically separate from finite-size operator identification.

Already established bounded results include:

- no exact finite nondegenerate **independent-bond** four-terminal gadget representation of one Bernoulli site in the declared class;
- no simple persistent bounded-degree factor mechanism in the explored exact matching polynomials;
- finite-size Galois complexity for several exact polynomials;
- the historical post-annihilator exponent near 7 is not automatically the next ordinary thermal spin-4 descendant.

Useful next work includes correlated-hyperedge/self-dual symbolic closure or obstruction and an honest post-annihilator correction-spectrum test. Both can run while Q1/Q2 work continues.

## Existing machinery that should be reused

### Acquisition optimizer — #102

This is **already implemented**. It successfully selected norm-5 as the high-information design for the previous H4/H12 model set. Do not list “build the optimizer” as future work.

Update its model/utility table to include:

- representation/projector separation;
- modulus/shape fingerprints;
- rank-revealing measurements;
- positive-control value;
- independence from already-used raw blocks.

### Sequential stopping — #126

Calibration already showed that the old H4-vs-H12 decision could often finish far below a fixed 100-batch budget, whereas H4-vs-zero was intrinsically much weaker. Future experiment design should therefore optimize **question + target + stopping rule**, not target geometry alone.

## Completed high-information blocks

- **#50 N145->290 full curve:** complete. Corrected slope/root structure survives; a single three-level multiplier shape does not.
- **#57 norm-5 N325/425:** complete. Frozen H4 beats H12/H8; child block alone remains compatible with zero.
- **#212 independent matching-odd synthesis:** complete. Global zero strongly disfavored; fixed H4 compatible.
- **#225 current multiradius prototype:** complete negative decision for the simple one-coordinate shell law; do not reinterpret this as a basis-invariant rank-2 falsification.
- **#248 deck-character selection rule:** exact. Global deck-invariant observables have zero first-order response to nontrivial deck characters; use charged/marked/twisted or nonlinear channels instead.
- **#260/#264 continuum spin-4 preflights:** complete as continuum arithmetic/velocity oracles; lattice field overlap remains the next question.
- **primitive square-bond sign/phase program:** complete for repeated norm-2 sign transfer; the next value is KdV/shape/character structure, not another same-purpose sign generation.

## Lifecycle debt to fix without rewriting history

The repository contains status drift between different views. A concrete example is the N145->290 block: current claim documents treat it as completed, while the primary evidence ledger still contains a `PENDING_REVEAL` lifecycle entry.

Do not edit a frozen preregistration simply to make dashboards agree. Regenerate/reconcile the lifecycle view from committed score artifacts and keep historical files immutable.

This is navigation maintenance, not a science gate.

## Decision logic

Choose the next experiment by the ambiguity it can kill and by actual maturity:

```text
lower spin-4 allowed or forbidden?       -> projector / charged selection work
known Jordan diagnostic works at all?    -> triangular positive control
H4 vs quotient/conjugation ambiguity     -> same-N coalescence #205
how many predictive states are needed?   -> #249 existing-data realization
specific Q4/Jordan shape vs alternatives -> typed modular-scalar modulus test
q2 vs Jordan dyadic composition          -> norm-4 #154 if acquisition ranking supports cost
why is p_c that microscopic number?      -> correlated-hyperedge / exact threshold-origin line
```

A failed discriminator is a successful result if it removes a mechanism class. A low-priority task is still allowed to run whenever it becomes cheap or yields a new information axis.
