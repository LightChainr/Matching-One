# Next Targets: Mechanism-Changing Decisions

**Updated:** 2026-08-29

This is the fast decision board. It asks:

> Which next output would most change the scientific model space, given the state, source, observer, geometry and acquisition already available?

It is not a queue or permission system. The seven lanes may run in parallel; priority changes attention, not authorization. Lower-ranked or unlisted work is not rejected, locked or required to merge.

## Read lifecycle before priority

One line may have a main theorem, a branch-only consequence, an open-PR result and a proposed production. Use these labels:

```text
main_integrated  canonical shared result
open_pr         inspectable but not merged
branch_only     inspectable remote frontier
hypothesis      proposed mechanism/discriminator
```

An old issue title does not override a completed result. A post-reveal fit is not a new independent block. Several derived coordinates from one histogram remain one dependency group.

## Seven current decision targets

| Rank | Decision target | Five-axis vector | Actual maturity | Exact next output | Acquisition |
|---:|---|---|---|---|---|
| 1 | **Why does the second activation change from reinforcement to cancellation?** | finite topology / Bernoulli thermal / global ordinary / Gaussian geometry and sector / reuse | The Draft K1/K2 reanalysis is complete: K1 dominates every size; K2 reinforces seven and has unresolved cancelling point estimates on three | Geometry/sector classification of the K2 component with existing covariance, followed by one held-out geometry only if required | existing data first |
| 2 | **Which complex C3 character is present?** | continuum/vacuum shape / typed response / complex topology row / three `rho` children / production | Exact child-character algebra is `branch_only`; current real rows alias conjugate characters | Three-point complex DFT with pure r=0/1/2 scores: scalar/constant, E4 weight-4, E4-squared weight-8, plus mixture residual | one typed complex row over three registered geometries |
| 3 | **One generator or context/morphism-enriched state?** | low-dimensional dynamics / common source / common readout / Gaussian×annulus rectangle / reuse+targeted acquisition | Norm-4 and annulus results exist off-main; shared rectangular realization is not yet built | Held-out covariance score comparing one shared generator with the smallest context-enriched alternative | reuse first; acquire only missing rectangle cells |
| 4 | **Does the log-pair statistic survive gauge and microscopic change?** | continuum Jordan control / cutoff shear / top-plus-two-spin and complex H4/E6 / triangular+square / targeted production | Triangular cross-cutoff result is `open_pr`; normalization audit and invariant design are `branch_only` | Gauge-invariant `CL` plus shear gate; then weight-12 two-modulus H4/E6 null in a common complex frame | small new control statistics |
| 5 | **Can the full Q derivative be reproduced on a known control?** | confluent module / Q measure+projector+field / boundary or tiny VJS / cross-ratio / exact+pilot | Exact measure, projector, Ward and boundary pieces exist across branch-only sources; no end-to-end calibrated vector | Anchored function-valued tangent reproducing the inhomogeneous ODE, `sqrt(3)/pi` log coefficient and required prefactor/field terms | exact construction, then smallest pilot |
| 6 | **Does the charged norm-5 sector have one fusion eigenphase?** | charged module / Z5 deck source / cubic charged row / N325 handed cover / engine reuse | Linear charged score is `branch_only` and H8-leaning; two primitive cubic channels are frozen in theory | Joint complex closure `C_A,+C_B,--C_A,-C_B,+=0`, normalization-free conjugate-paired magnitudes and full covariance | reuse engine; new marked three-point accumulators |
| 7 | **Where can a nontrivial relative-source deformation live?** | defect/connectivity radical / relative fugacity / junction observer / exact finite algebra / theory | Scalar source closes at rank three and its semisimple algebra has `HH^2=0` on branch-only proofs | Minimal enlarged connectivity algebra, separable quotient, radical/module-extension data and one explicit nontrivial or obstructed class | exact work only |

## Target 1 contract: two-activation H4 decomposition

The exact definitions are

```text
K1 = K_minus,
K2 = K_plus,
F1(p) = E[H_K1(p)],
F2(p) = E[H_K2(p)],
M(p)  = -1 + F1(p) + F2(p),
C = (K1+K2)/2,
G = K2-K1.
```

For each aligned direction pair at common root `p_bar`, compute

```text
DeltaF1, DeltaF2, DeltaM=DeltaF1+DeltaF2,
delta_p1=-DeltaF1/Mbar',
delta_p2=-DeltaF2/Mbar',
actual root difference - (delta_p1+delta_p2).
```

Use only exact `DeltaCos4` direction normalization; do not fit a free exponent. Delete one aligned batch across both directions, not one component at a time. Shared streams remain one dependency group.

Default archives are N=65,85,130,145,170,185,265,290,325,425. A dataset without the required `K_minus/K_plus` semantics or aligned batch structure must be reported as `not_scoreable`, with a reason; do not substitute a similar-looking histogram.

Entry point on this Draft branch:

```bash
python3 scripts/analyze_two_activation_h4.py \
  --manifest analysis/two_activation_h4_manifest.yaml \
  --output-json results/two-activation-h4/latest.json \
  --output-md results/two-activation-h4/latest.md \
  --workers 16
```

Completed result:

```text
scoreable sizes: 10/10
larger component: K1 at every size
K2 reinforcing:   N=65,85,130,145,170,185,290
K2 cancelling:    N=265,325,425 (each negative term individually |z|<2)
```

The nonlinear closure residual is negligible relative to the root gap at every size. The next question is therefore not whether the linear decomposition works, but which geometry/sector feature controls the sign and relative size of `K2`. N=65,85,130,170 remain one shared dependency group.

Interpretation outcomes:

| Pattern | Scientific update |
|---|---|
| `DeltaF1` dominates | H4 is primarily first appearance of nontrivial homology; this is the observed point-estimate pattern at all ten sizes |
| `DeltaF2` dominates | H4 is primarily completion of the second ambient direction |
| large opposite components, small sum | the matching root hides a topology-sector cancellation |
| large closure residual | linear root-shift decomposition is insufficient at current anisotropy; retain the exact curve-level split |

## Target 2 contract: complex C3 character

Use one typed complex/chiral observer and the three degree-2 child geometries in one registered frame. The frozen ordinary-ring predictions are

```text
r=0: (1,1,1)             E6 / constant child character
r=1: (1,zeta,zeta^2)     E4 / weight-4
r=2: (1,zeta^2,zeta)     E4^2 / weight-8
```

Keep real and imaginary covariance together. A real-only `1:-1/2:-1/2` score cannot distinguish r=1 from r=2. More than one DFT component is an interpretable mixture/tangent/nonlocal result, not an automatic failure.

## Target 3 contract: context rectangle

The comparison is not “fit another rank everywhere.” It must use the same declared source/readout basis on a rectangular subset of:

```text
contexts = {Gaussian scale/Smith/deck, annulus radius/geometry},
outputs  = {ordinary or typed matching row, shared covariance coordinates}.
```

Fit on a strict subset and hold out at least one Gaussian and one annulus cell. Compare:

1. one common generator/state realization;
2. the smallest context- or morphism-enriched realization that adds an intermediate filtration/label.

The final CRT endpoint cannot be used as a path-order label. If chronological memory is claimed, retain intermediate ranks or dynamic lineage.

## Target 4 contract: gauge closure

Do not treat the triangular `kappa_proxy` or a bare square-H4/triangular-E6 amplitude ratio as universal. The next rows must cancel the relevant gauges:

- normalized top-plus-two-spin statistic for `CL` and its shear ratio;
- two-modulus, weight-12 H4/E6 polynomial null across the same C3 children;
- optional `LL,LD,DD` Gram invariant only with its finite-cutoff shear boundary stated.

The current open-PR result already establishes a nonzero cutoff shear within its declared gauge. Repeating the same two-point production without a new normalization axis is lower information.

## Target 5 contract: Q-tangent control

The saved comparison vector must separately expose:

```text
measure score,
finite confluent projector contribution,
explicit field/insertion derivative,
boundary conformal-prefactor derivative where applicable.
```

Use the Ward `4:-6:3` row to test thermal descendant versus primary and the projector residue direction to test singlet/`[2]` confluence. For the boundary route, score the anchored cross-ratio tangent after removing the single amplitude gauge. Passing only the measure covariance is not an end-to-end control.

## Target 6 contract: charged fusion

Charge conservation leaves the primitive unordered triples `113` and `122` and their conjugates. Their one- and two-point disconnected pieces vanish exactly. Preserve the transported deck basis and save the full covariance of the four complex handed channels.

The existing `71/21/8` ordering is one 4D comparison, not four votes and not a spin-8 discovery. The cubic continuation asks a new operator-algebra question and remains distinct from the ordinary global H4 block.

## Target 7 contract: radical rather than another scalar source

Do not add `q^3,q^4,...` as nominally new source coordinates: they reduce exactly to the three rows `1,q,q^2`. Do not seek a Hochschild class by merely changing coefficients or attaching an ordinary bimodule to the same separable scalar algebra.

The next algebra must retain information absent from the scalar sectors, such as connectivity partitions, junction composition, non-semisimple radical states, doubled-space exchange, singular projectors or marked/charged defect endpoints.

## Secondary opportunity pool

These lines remain active and may overtake the default order whenever they obtain a sharper discriminator or low-cost acquisition:

| Line | Current boundary | High-value continuation |
|---|---|---|
| Norm-4 / norm-10 | q2 rejected; Jordan plus one even mode is post-reveal economical fit | freeze the source-only matrix before a genuinely held-out Gaussian target |
| same-N norm-5 coalescence | exact amplitude-free interpolation is ready | use when its new Smith class fills the context rectangle efficiently |
| local pivotal | current N130/N170 rows nearly rank one | add a thermal-null or representation-typed row, not merely replicas |
| modular-scalar modulus | channel typing is `main_integrated` | one new typed `cross/either` shape with a parameter-free prediction |
| microscopic threshold origin | exact semantics/evaluators exist; no threshold bound | correlated-hyperedge, block-event confidence or symbolic obstruction |
| post-annihilator corrections | composite algebra rules out some shortcuts but not amplitudes | joint harmonic/phase sidebands rather than exponent alone |

## Recovered frontier sources

| Node | Lifecycle pointer |
|---|---|
| digital Alexander proof | `main_integrated` `2da5855`, `notes/digital-alexander-duality-proof.md` |
| two activation | `branch_only` `theory/p28-two-activation-rank-mixture-20260829@b8004bc`, `notes/two-activation-rank-mixture.md` |
| global selection | `branch_only` `analysis/p257-global-singlet-selection-20260829@9320649` |
| small-width selection | `branch_only` `analysis/p120-operator-spectroscopy-20260829@9cf2373` |
| norm-4 reveal / fit | `open_pr` #273 `8b26a30`; #277 `039e708` |
| annulus held-out recurrence | `branch_only` `analysis/p253-annulus-recurrence-20260829@8e91c90` |
| ordered filtration | `branch_only` `analysis/p200-path-ordered-filtration-20260829@fb82de4` |
| triangular log pair | `open_pr` #246 `7f9dcd8`; normalization audit `branch_only` `ceb7c6e` |
| charged norm-5 | `branch_only` `analysis/norm5-chiral-hecke-phase-20260829@cc1d43c` |
| Q/projector/Ward controls | `branch_only` `1055e22`, `d006f9c`, `fa73f5d`, `2a15e63` |
| C3 character and invariants | `branch_only` `57b59be`, `a6120aa` |
| relative-source closure/rigidity | `branch_only` `94dd7f9`, `094ee36` |

See `docs/STATUS.md` for full branch names, artifact paths and claim boundaries.

## Stop/re-rank rules

- Do not count a post-reveal coordinate on an old stream as new primary evidence.
- Do not increase N365 alone when uncertainty is dominated by old recurrence coefficients.
- Do not add replicas to the same N130/N170 local pivotal rows when a new observer is the missing dimension.
- Do not infer path order from a final unmarked CRT join.
- Do not pool ordinary and charged norm-5 rows as the same operator amplitude.
- Do not call a gauge-fixed agreement universal before a normalization-free statistic exists.
- Do not use the roadmap to forbid a mathematically distinct attempt; re-rank it by the new information axis it supplies.

The target is a faster sequence of scientific separations, not a smaller number of ideas.
