# Research outputs: thermal windows, named readouts, and a physical cut state

This update follows the [#469–490 context increment](../docs/ISSUE-PR-INCREMENT-20260831b.md)
and an additional bounded read of PR #491. It does not restart the full-repository
review or the completed tests. Two new analyses consume existing data below.

## N100: the center is not silent, even though its odd dipole is unresolved

Open [PR #484](https://github.com/LightChainr/Matching-One/pull/484), pinned at
`894b3d800c5aeaad3dd8b0f893b6f17d85d234c6`, has already produced the three
N100 shapes and rejected a common affine scalar profile for A/E/C/W. Its later
finite-Jacobian transport comparison is also complete: six cumulative
invariants give nominal chi-square 53.91436/6. A alone admits an empirical
quantile map; the shared A/E map is the discriminating question. Neither the
first production nor another low-degree warp is still pending.

Draft #267 commit `d973a39` now locates the clock-quotient residual
`R=U-r_C D` in fixed regions, using the exact saved marginal histograms.
For `z=N^(3/8)(p-p_ref)` the center `|z|<=1` is
`p in [0.414918110,0.770573992]`. This is a declared finite-N coordinate,
not a fitted exponent or a proven asymptotic critical window.
[Report](../results/etop-critical-window-n100/REPORT.md),
[complete covariance and all 200 deletions](../results/etop-critical-window-n100/latest.json).

| Center-window readout | Integral against dp | Aligned delete-one SE |
|---|---:|---:|
| R_A | −3.460293249e−4 | 4.626773252e−5 |
| z R_A | −1.011331359e−5 | 2.241123543e−5 |
| z² R_A | −8.532091396e−5 | 5.559638025e−6 |
| R_E | 1.619451641e−6 | 5.549053422e−5 |
| z R_E | −1.358016369e−4 | 1.167770404e−5 |
| z² R_E | −7.650714175e−6 | 1.528835254e−5 |

The full odd dipole is −2.93635e−4. Its lower/core/upper pieces are
−4.19879e−4, −1.01133e−5 and +1.36357e−4, with their joint covariance retained.
These are **signed additive moments**, not positive mass or causal shares.
The central A area and E dipole are resolved while their complementary moments
are not. Thus “the full dipole is mainly outside this window” must not become
“the center has no response.” Matching parity (A/E) and parity about p_ref
are also different; this finite pattern is not an exact thermal-parity law.

In first/second-activation coordinates, the central areas are
`R_F1=-1.738243883e-4 +/-3.50209e-5` and
`R_F2=-1.722049366e-4 +/-3.71953e-5`.
They reinforce the central A area; their difference gives the unresolved E
area. This is a linear change of coordinates on the same random block,
not two new replications or an operator assignment.

All three shapes use the same 2M permutation counters and 200 batches from
PR #484. PR #485's clock quotient, the full-p transport tests and this window
analysis are dependent. Every deletion removes the same batch across all
shapes and orientations and refits r_C. There is no new MC; incomplete-beta
identities are evaluated numerically, not by p-grid quadrature or a certified
interval method. Windows .5/1.5 and all birth moments are included as correlated
retrospective descriptions, not additional significance votes.

**Next useful prediction:** carry the same moment vector across homothetic
scales, or first allow a separately estimated common center-window secant to
distinguish global-clock calibration from a truly extra central response.
PR #484 already contains a future N400 design; its manifest is not acquired
N400 data. The original norm-4 local energy/singlet identification remains a
distinct question and cannot be declared complete from these topology curves.

## P398: an equal-time-defined observer exposes a real fast response

Draft #267 commit `4846adf` consumes the existing width-five C(d) and residues.
The landing innovation

`J=L-[C_LA(0)/C_AA(0)] A`

is chosen using equal-time covariance only, not a nonzero-lag fit or a spectral
null. At d=1 the unit-variance signal is **.0202159964** and the fixed two-slow-mode
approximation misses **77.9824%**, compared with .07698% for A, 1.97589% for L,
and .8996% for the whole matrix. [Report](../results/p398-fixed-readout/REPORT.md).
At d=4 the large relative J error accompanies an absolute signal of only
3.30e−8; quoting the percentage alone would be misleading.

The additional source/future mode-null filters are explicitly model-derived
exploration. Their near-100% two-mode loss is constructed, not a separate
falsification. The equal-time innovation is the useful independent selection
rule within this same finite model. Experimental detectability would still
need sampling and coefficient-calibration covariance, which this exact-model
archive does not supply.

New branch-only
[`552c45d`](https://github.com/LightChainr/Matching-One/blob/552c45d7595ebcb0d04555cec03b2a5bfd8da44a/notes/p398-width8-source-spectrum.md)
has separately completed width-eight continuous propagation with fixed i
character. Kreweras symmetry protects two readout rays, but each generates
93 directions. The microscopic leakage is explicit:
`GL=-3L+T2`, `GA=-3A+R`, with size-two membership T2 and boundary-contact
multiplicity R. These are not the width-five discrete ζ5 filters.

The next step starts **after** first fixed-readout/width-eight propagation:
use those geometric emissions and protected rays for a common-model width or
microstructure prediction. Exact state dimension, useful compression for a
named observer, and continuum field identification remain separate objects.

## P334: the two graph sides now have an occupied-cut construction

Open [PR #491](https://github.com/LightChainr/Matching-One/pull/491), pinned at
`ab90201e88409310632812727e0138c56b455644`, supplies the
[proof and saved-checkpoint explanation](https://github.com/LightChainr/Matching-One/blob/ab90201e88409310632812727e0138c56b455644/notes/p487-cut-network-theorem.md).
For the stated embedded rank-one torus scope, cutting an occupied essential
cycle turns second-rank continuation into two-terminal **vertex** connectivity.
The sides arise before measuring pair edges. Neutral occupied components
generate bicliques, while direct vacancy edges supply the other possible term.

The two old N425 witnesses are explained as K(12,8)+K(3,4) and
K(25,4)+K(4,2), each with one shared vertex. Their W2 difference is
`540=472+68`, and the same contact networks give c3=583/509.
The fixed-cut network also closes under activation/contraction. This is not
a proof for crossed matching graphs, a minimal bounded state, a population
H4 contribution or a new independent sample block.

The useful continuation is now a cut-invariant/covariant component-incidence
response on the existing paired population, or a controlled longer-horizon
vertex-reliability prediction. Do not repeat the first cut implementation,
the two graph reconstructions or the already explained W2 collision.

## P418: correct the archived comparison, then change the scientific interpretation

The same radius4/5/6 production CSVs store batch sums with exposures
200/3000/3000. The original common reader did not divide by exposure.
Draft #267 `e2b57aa7`'s [normalized archive result](../results/p418-normalized-archive/REPORT.md)
changes only that input unit before the inherited geometry/mask, covariance
and family fit. No old input or score is overwritten.

The corrected common masked distances are about70.878/75.692/76.766/68.067,
with nominal bootstrap p-values .382/.183/.175/.367 under the inherited
250-draw convention. The saved nonnegative solutions attain the resolved-rank
linear least-squares lower bound to below7e−11. The historical584–1153 common
masked distances therefore do not support the previous radius-flow story.

There is a separate numerical limitation: radius5 alone has no zero-lag
mass anchor, and the inherited rank-deficient NNLS uses enormous near-null
weights. Its separate-shell distances and the sharing penalties that subtract
them are marked numerically unreliable, not repackaged as a new mechanism
result. The common fit includes radius4 zero lag and passes the point-solution
residual/optimality diagnosis. No second bootstrap or full test suite was run.

Common raw and masked compatibility is not uniqueness of a reconstructed
spectrum or proof that the two cones are universally equivalent. The exact
CRT/root-translation results, independent normalized P250 state/rank analyses
and paired-anchor pilot's own statistics are retained. A genuine separating
spectral/readout prediction is now more informative than explaining the old
unit-induced large penalty with a new physical mechanism.

## Parallel scientific work

Use the [single attention board](../docs/NEXT-TARGETS.md). The new six #370
framework PRs and eight finite serial-monoid results are available support,
not reasons to return all teams to preparation. One team can pursue the
original norm-4 independent physical source, another the completed thermal
profile's scale prediction, and another the now-explicit connectivity state.
These are parallel suggestions, not locks or approvals.

All new local calculations above reused saved inputs. No server/tunnel,
new production, test suite, Issue close/lock/rename, merge, rebase or force-push
was performed in this continuation. New branches are read and cited with
their integration status; they are not silently merged into main.
