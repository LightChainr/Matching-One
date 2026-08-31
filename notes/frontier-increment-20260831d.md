# Frontier increment d: completed physical analyses after the context review

**2026-08-31. Read and organize; no new calculation.** The substantial change is
not another preparation layer: the full continuation calculation now covers
147 fixed production prefixes; its canonical ordering and removable suffix
noise have actual results; a stationary-current deletion has been evaluated;
and the thermal-shape analysis excludes a broad symmetric two-lobe mechanism.

This note follows [context reconciliation c](context-reconciliation-20260831c.md)
and feeds the existing [Next Targets](../docs/NEXT-TARGETS.md). It does not
introduce a second queue, stop parallel work, or make tool completion a
condition for research.

## Read boundary and integration

The bounded GitHub read at **2026-08-31 06:41:26 UTC / 14:41:26 Asia/Shanghai**
captured open, non-Draft PR #484 at
`8b3c4e4b175826805238f2585ca470d746a80610`, and
[PR267 comment5474724807](https://github.com/LightChainr/Matching-One/pull/267#issuecomment-5474724807).
The comment located the cross-branch results; it is **not** the basis for
claiming their full review. The four requested P334/P398 scientific notes and
their reports/cards were then read completely from their existing local Git
objects. The two new #484 notes and thermal-model report were also read
completely. No network solver, Monte Carlo, scientific test or result script
was executed.

| Scientific result | Integration / branch | Reviewed immutable commit |
|---|---|---|
| Common symmetric two-lobe obstruction | `open_pr #484`; `analysis/etop-modulus-survivors-20260831` | `ddf7d564cb2a99bce4bbdba5eef90a9143e7d394` |
| Log-occupation race and initial acceleration | `open_pr #484`; same branch, captured head below | `0997a00c4639af90b806b6920489e32cc1e85227` → `8b3c4e4b175826805238f2585ca470d746a80610` |
| All147 full physical clocks | `branch_only`; `analysis/p334-all147-prefix-clocks-20260831` | `87b6ca5b39084c06143f31cafdaba53f90012e27` |
| Fixed147-mixture noise budget | `branch_only`; `analysis/p334-147-prefix-noise-mixture-20260831` | `5a8e26ec8cd63fb786e212db601c242d882c92bc` |
| Twelve-prefix canonical crossings | `branch_only`; `analysis/p334-twelve-canonical-crossings-20260831` | `1b4f5499b68c0b2a3881108b2219f703128d6594` |
| Width8 stationary-current deletion | `branch_only`; `theory/p398-width8-reversible-current-control-20260831` | `520a9d218a90ce46c34475e2ddd4ef0eb5d97e1a` |

The three P334 cached remote refs matched their reviewed commits. The P398
cached ref had already advanced to
`9339fe7053e7edb9573431eb8dc254be91bfb388`; that later content was **not**
followed in this bounded read. Its results must not be inferred from the ref
name or folded into the reviewed520a9d21 result. The table records scientific
pins, not a claim that each pin remains the latest branch head.

Cross-branch citation in #267/#484 does not make those files integrated in
either PR. In particular, the standardized-rank-shape input
`6d8a3ed9d961c66889c3c1e4575485443fdd1c39` belongs to the separate
`analysis/p267-scalar-clock-transport-20260831` branch. Its consumption by
ddf7d564 does not silently import that source branch or create new data.

## 1. Thermal shape: varying two peak positions and weights is insufficient

Read: [mechanism note](https://github.com/LightChainr/Matching-One/blob/ddf7d564cb2a99bce4bbdba5eef90a9143e7d394/notes/etop-two-lobe-moment-obstruction.md),
`results/etop-two-lobe-moment-closure/REPORT.md` and `score.json`.
The companion is `scripts/etop_two_lobe_moment_closure.py`; it was inspected,
not run.

The observer is the full-interval signed rank-step difference
`D_A=Y_(4i)-Y_(2i)`, standardized by its own moments. A positive two-component
mixture is a **candidate representation**, not an assumption that this signed
response must already be a probability density.

For equal-width Gaussian lobes, moments3/4 determine center separation and
weight; moments5/6 give nominal joint residual848.387/2 at N100 and177.594/2
at N400. The stronger result removes the Gaussian assumption. For
`X=B+Z`, with B taking two centers and Z any common independent symmetric
positive kernel, moments3/5 determine the centers. The remaining kernel
moments must satisfy positivity:

| Necessary quantity | N100 | N400 |
|---|---:|---:|
| Kernel sixth moment | −1.757800 ± .074189 | −3.849346 ±1.20234 |
| `v E[Z^6]−E[Z^4]^2` | −.161132 ±.008132 | −.297007 ±.050780 |

Thus the current mean moment vectors cannot be represented by two translated
copies of one symmetric positive kernel, even allowing a different kernel
at each N. Unequal peak shapes, within-peak asymmetry, more components or
signed cancellation remain distinct possibilities. This is not merely a
failure of one Gaussian fit.

Every source-aligned delete-one vector refits the implied centers and
retains joint parameter/residual covariance. N1002M and N4008M are independent
source blocks; Gaussian and non-Gaussian scores within a size are dependent
views. The model was chosen after source reveal. N400's inferred kernel
variance is near zero, and propagated errors are not a calibrated boundary
likelihood test or a strict finite-sample confidence certificate.

**Completed, not next:** first two-Gaussian closure and first arbitrary-common-
symmetric-kernel obstruction. The frozen N900 width question remains different
from this later exploratory full-shape question.

## 2. Full birth law: independent competition has an exact clock

Read the complete
[general two-port note at8b3c4e4b](https://github.com/LightChainr/Matching-One/blob/8b3c4e4b175826805238f2585ca470d746a80610/notes/p334-general-two-port-birth-theorem.md).
Its two-port theorem and safety-polynomial factorization were already in the
previous read; the new contributions here are the logarithmic-time race and
short-time acceleration.

For d vacancies and h original singleton triggers, transform independent
uniform occupation marks by `s=−log(1−u)`. The direct clock is Exp(h) and is
independent of the collective completion event on the other sites. Therefore

\[
\pi_{\rm direct}
=h\,E[s_{\rm birth}]
=h\,E[H_d-H_{d-T}],
\qquad
\operatorname{hazard}_{\rm birth}(s)
=h+\operatorname{hazard}_{\rm collective}(s).
\]

The original-direct fraction is a harmonic moment of the **complete** discrete
clock, not a function of H2 alone and not a new sitewise-solve requirement.
Independence holds for this transformed independent-mark construction, not
for component clocks conditioned on a fixed total insertion count.

With m2 safe-site minimal pairs, the initial log-time hazard derivative is
exactly2m2. For counter suffixes83/1006,
`log(S83/S1006)=−s+18s²+O(s³)`. The resulting short-time estimates4.74/9.35
expected insertions explain the separation between the observed discrete
hazard reversal at5/6 and survival reversal at10/11. They are not finite-step
error bounds.

The topology remains the genuinely embedded black-NN rank-one torus graph.
It does not cover unjoined crossing diagonals or imply generally small
network treewidth.

## 3. All147 clocks: irreducible triples break a stronger scalar closure

Read:
[full note](https://github.com/LightChainr/Matching-One/blob/87b6ca5b39084c06143f31cafdaba53f90012e27/notes/p334-all147-real-prefix-clocks.md)
and `results/p334-all147-prefix-clocks/SCIENTIFIC_CARD.md`.
The result directory preserves `selection.json`, `full_clocks.json`,
`prefixes/<counter>.json` and `scientific_summary.json`.

The source selection was frozen at
`8d7ac0e91323608053df3a8d8e06db02335fb315`: all147 eligible rows of the
existing N425 second-orientation source at k0=252, age10, ell=(12,−19),
excluding the original two case-study counters. The old12 are reused;135
additional conditional networks are evaluated. The saved
`9cca7bc60e26db5ec47b5e00fbc5d98532447c29:results/p334-all147-prefix-clocks/full_clocks.json`
records `new_samples=0`,135 new rows,12 reused rows and14.45890854 seconds
batch wall time. That is a stored execution record, not a runtime estimate
made by this review.

All147 solve under the fixed caps, with535 factors, maximum factor treewidth9
and maximum boundary-state count25,882. Mean waits range5.77269–30.86927.
Of10,731 pairs,3,510 have exact survival-order crossings;1,867 first switch
by step40. These are properties of the fixed eligible source stratum, not a
population crossing frequency or universal tractability bound.

The strongest mechanism witness is more specific than a large clock spread:

| Same N, geometry, age and line | Counter43042508631 | Counter43042514803 |
|---|---:|---:|
| H2 / safe pairs / safe-successor-square sum |15 /12,397 /3,890,796|15 /12,397 /3,890,796|
| True safe triples |644,020|644,006|
| Genuine minimal triple triggers |5|19|
| Mean birth step |10.48090361|10.39836128|

The stored second moment fixes minimal-pair count6 and wedge count5.
Together with the embedded black-NN bipartite-trigger result, it yields
pair-only safe triples644,025; subtracting the actual coefficient exposes
irreducible triples5/19. No fresh triple census is needed. Their first
difference already occurs at insertion3: second-moment cooperative branching
does not close this clock. This does not attribute the whole late-time gap
only to triples or establish path memory.

The first complete-prefix expansion, fixed-stratum cost observation and
direct/collective shares are **done**. The scientific card's proposed first
thermal/noise readout is superseded by the following completed results.

## 4. Canonicalization retains the twelve states' ordering reversals

Read:
[canonical note](https://github.com/LightChainr/Matching-One/blob/1b4f5499b68c0b2a3881108b2219f703128d6594/notes/p334-twelve-canonical-order-crossings.md)
and `results/p334-twelve-canonical-crossings/REPORT.md`; `score.json`
contains66 pair certificates and21 order cells.

All11 crossing pairs among the original fixed12 retain their rank reversals
after the canonical binomial readout, producing20 simple thermal roots;
the other55 pairs remain ordered. One coefficient sign change forces exactly
one simple interior root by Bernstein/Descartes variation. The survival of
both roots for the nine two-change pairs is an actual finite certificate,
not a general theorem for all two-change profiles.

For `Delta=F2_83−F2_1006`, the unique crossing is p=.594353897611717.
At p_ref=.59274605079, Delta=+.0002979084441, while its integral is
−.00275241089307. Thus reference-p ordering and integrated waiting-time
ordering can disagree on these same conditional states. Canonicalization is
not the substitution p=(k0+k_cross)/N.

The root counts are certified by exact reduced Bernstein coefficients and
dyadic partitions; locations/slopes are numerical evaluations. Mathematical
roots need not be useful measurements: the48/622 root near.996 has a remaining
tail around6.71e−46. All p views reuse the same12 conditional laws, which are
included in the147 set; no independent samples are created.

## 5. The first147-prefix noise budget is already measured

Read:
[noise note](https://github.com/LightChainr/Matching-One/blob/5a8e26ec8cd63fb786e212db601c242d882c92bc/notes/p334-147-prefix-conditional-noise-weight.md)
and `results/p334-147-prefix-noise-mixture/REPORT.md`; `score.json`
retains full two-readout within/between/total covariance and all147 rows.

The defined experiment chooses one of the147 prefixes **uniformly**, then
a fresh uniform suffix. Its readout is the real binomial tail
`g_T(p)=Pr{Bin(425,p)≥252+T}`, not a Bernoulli draw, and
`integral g_T=(174−T)/426`.

| Readout | Mixture mean | Suffix variance | Between-prefix variance | Removable fraction |
|---|---:|---:|---:|---:|
| g_T(p_ref) |.1706204126|.02092217261|.003543006813|85.5182%|
| Integrated g_T |.3744119018|.0005775498650|.0001104170921|83.9502%|

This is the law-of-total-covariance decomposition of the declared finite
mixture, with denominator147, not a sampling ddof correction. Integrated
quantities also have exact rational storage. If all147 conditional means
are averaged directly, that finite-mixture mean is already deterministic.

The result is a concrete existing-data opportunity for conditional averaging,
but **not** a global Monte Carlo speedup, stratum-frequency estimate,
cross-orientation covariance or full A_top/E_top result. K1 is not reconstructed
from age. Prefix generation cost and coverage beyond this stratum remain
different questions. No new network solve, random suffix or p-grid was needed.

Do not reassign the first mixture variance calculation after reading the
older147-clock or twelve-crossing note's “next” paragraph.

## 6. P398: irreversible current amplifies the tail, but does not create inversion

Read:
[current-deletion note](https://github.com/LightChainr/Matching-One/blob/520a9d218a90ce46c34475e2ddd4ef0eb5d97e1a/notes/p398-width8-current-deletion-result.md)
and `results/p398-width8-reversible-current-control/SCIENTIFIC_CARD.md`.
The fixed-grid numerical outputs, spectra, residues, hashes and arithmetic
identities are in `latest.json`; the definition was frozen at
`5e47bdb643a608d8086bd30a28b6f7044d4ea5b1`.

The single counterfactual replaces G by S=(G+G*)/2 on the same1430-state
positive width8 connectivity process. It preserves stationary pi, fixed
AP/landing rays, source variances, every state's exit rate and real initial
slopes. It removes every stationary probability current.

| Fixed observer quantity | Original G | Current-deleted S |
|---|---:|---:|
| Normalized ray crossing |.2656573200|.2722634760|
| Lowest visible mass, minus |2.8196586326|2.5407959794|
| Lowest visible mass, plus |1.9557501384|1.8363180504|
| Plus/minus correlation ratio at s=4 |25.41893805|14.42424134|

The plus ray still starts faster and ends slower. Ordinary reversible
positive mixtures suffice for this inversion; negative residues, complex
modes and a Jordan block are not necessary in this finite example. Current
amplifies the long-time contrast while moving the crossing only2.4867%.

Deleting current increases initial log curvature by
`||J psi||_pi²/||psi||_pi²`, .3808266892 for minus and.3600854185 for plus.
Their similar absolute contributions explain why the early crossing changes
little despite appreciable individual memory changes. G and G* have identical
real self-correlations, so those correlations cannot determine circulation
direction. The520a9d21 note identifies an allowed anti-Hermitian
cross-correlation as the directional object; it does not compute that later
readout. No claim is made here about the unreviewed9339fe70 continuation.

This is a completed deterministic finite-process intervention on the same
source as the width8/memory/motif results, not an independent data block.
S is a legitimate reversible generator but need not be a local square-bond
transfer word; it does not identify a site-Matching field or universality class.

## N900: distinguish authored progress from a completed artifact

At the captured #484 head, the new body reports that two local single-thread
N900 processes are running, with32M shared counters/800 batches and two shapes.
It expressly reports **no target result yet**. Compared with705819e, the
three new commits add only the two-lobe result and two-port-note extensions;
there is no N900 raw/result/completion receipt in that increment.

Accordingly this note records **author-reported running; completion not
established; process state not independently inspected**. The earlier c note's
artifact-only/runtime-unknown entry remains its historical observation. The
existing frozen width predictions2.565535/2.094751 share the N400 anchor
covariance; target uncertainty is additional. This is not authorization to
restart, duplicate or stop production, and no process/server action occurred.

## Consequence for the single research handoff

The original norm-4 physical response remains the leading unresolved identity
question; none of the finite-clock or reversible-frontier results names that
field. They do remove several concrete preparatory interpretations:

- **Clock structure:** use the complete147 laws, marked competition and actual
  mixture noise budget. The missing scope is genuinely paired/population and
  stratum-weighted information, not a first extra-prefix, first full clock,
  first variance decomposition or another generic triple count.
- **Shape:** consume the existing N900 campaign when its result exists. Keep
  its frozen width decision separate from a post-source exploratory test of
  asymmetric/unequal/signed lobe structure; another weight/separation scan
  cannot repair the measured common-symmetric-kernel obstruction.
- **Propagation:** current deletion is completed, not the next first
  intervention. A specifically defined microscopic/matching overlap or new
  physical generator response remains distinct from this finite symmetrized
  control; later work on the same branch must be read before assigning it.

All P334 views retain the old N425 source dependency, with12 nested in147.
All P398 views retain the same finite generator/source dependency. The
N100/N400 shape views retain their two original blocks. No source is added to
an independent evidence total merely because a new report or coordinate exists.

Only this note was written. No Issue/PR mutation, merge, history rewrite,
simulation, scientific test or server operation accompanies it.
