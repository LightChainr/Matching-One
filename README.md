# Matching One

Matching One studies square-lattice site percolation through its matching-lattice identity, finite topology, operator sectors and the microscopic origin of the threshold. The exact structural anchor is

\[
p_c^{\mathrm{site}}(\mathbb Z^2)+p_c^{\mathrm{site}}(\mathrm{NN+NNN})=1.
\]

The repository is organized to expose the next mechanism-changing observation, not to turn priorities into permissions. Exact work, reanalysis, pilots and independent theory lines may proceed in parallel. A lower priority is not a rejection; no task is locked by the overview.

## Start here

**Current decision: #537 has crossed the finite transmission gate and moved
to P1 theory; #275 is the new P0 attention line.**  Finite pure-thermal rank
one and the automatic two-independent-defect/six-arm story are retired.  The
exact N25 survivor is a nonfactorizing `contact fusion x topological
completion` tensor, and a frozen 20M held-out N65 block reproduces its four
signs `-- -+`.  The same sufficient statistics now also give the complete
canonical response `J65=-.00162250989+/-.00018553008`.  The typed carrier
share falls from `5.892%` of full `T25` to `2.551%` of full `T65`: it is a real
transmitted microscopic coordinate, but not the dominant complete response.
Highest attention therefore returns to the original Matching-One observable:
state exactly what is traced, how it is normalized, and make at least two
candidate mechanisms predict distinguishable vectors on existing raw
coordinates.  No new Monte Carlo is the default.

The [frozen N100/N400 pilot](results/regular-pair-macro-joint-u/REPORT.md)
uses the same `Kreg`, the exact Euclidean window
`L/4 <= d_T(x,y) <= 2L/5`, paired occupation/anchor streams, exact
K-likelihood reweighting and the complete direct/root/source-slope/root-slope
original-U functional.  It gives

| extensive response | point estimate | simultaneous familywise-95% interval |
|---|---:|---:|
| `T100=100^2 J2_macro` | `+11.8778461` | `[-2.2556162,+26.0113084]` |
| `T400=400^2 J2_macro` | `-542.5038231` | `[-1307.642507,+222.634861]` |

Both intervals contain zero, the point estimates have opposite signs, and
the frozen N400 production projection is `2,119,100`, above the two-million
ceiling.  The declared stop rule therefore forbids a top-up or replacement
window and leaves `D17/D21` unevaluated.  No GPU or cloud job was needed:
100,000 new paired configurations took about 3.3 producer seconds in total.
The `s=2` coordinate carries `11.465/11.878` at N100 and
`-540.463/-542.504` at N400, but it shares every configuration and covariance
with `s>=3` and is not a second vote.

The route was not stopped because raw interaction vanished.  PR #509
[`eed2190c`](https://github.com/LightChainr/Matching-One/blob/eed2190c04b67084ab5aef5827e00377853a0bca/notes/p337-critical-spatial-summability.md)
proves at the exact square-site critical point that the raw canonical kernel
is absolutely summable, `E|g_xy| <= C d(x,y)^(-2-eta)`.  The
[thermal-pivotal analysis](notes/regular-pair-thermal-pivotal-gate.md) shows
why this does not close original U: an arbitrarily remote merger pivotal can
change `g_xy` from `1/16` to zero, and the direct/source-slope terms require
q/E-weighted three-point influences. Existing rigorous `alpha4>1` is enough
for raw summability but not for their position sum; the
[literature bridge](notes/p337-thermal-pivotal-literature-bridge.md) shows that
even the triangular-lattice value `alpha4=5/4` leaves the naive three-centre
absolute account growing like `R^(1/4+o(1))`.  The first finite rank test has
now split into a result and a correction.  PR #509
[`ec3941b0`](https://github.com/LightChainr/Matching-One/blob/ec3941b03b2694e827db1cba34766a82e6146a5a/experiments/p537-landing-matrix-preflight-20260901/REPORT.md)
finds all projected minors nonzero for an explicit provisional N25
six-block clean-two-bridge contract.  This kills that port-level rank-one
model.  The N9/N16 calculations at
[`b38fbb61`](https://github.com/LightChainr/Matching-One/blob/a3bc80c86585220f96dfeff022bd575b6c21d29f/notes/p537-ordinary-landing-minor.md)
and `2c4942b4` remain exact raw controls, but the corrected branch head
[`a3bc80c8`](https://github.com/LightChainr/Matching-One/blob/a3bc80c86585220f96dfeff022bd575b6c21d29f/notes/p537-finite-landing-transfer-definition.md)
shows why they do not finish the canonical gate: their columns split one
covariance identity, they do not contain separately normalized axis/tilted
P4, and the row shear omits the bilinear `-beta H B` allocation.

Closed-unmerged Draft
[#544](https://github.com/LightChainr/Matching-One/pull/544) preserves a
full-root/Schur **coarse-marginal** collar value test.  Its theorem at
[`139c7e58`](https://github.com/LightChainr/Matching-One/blob/139c7e5850c5be1e8e17c3d58be806080b1e2b73/notes/p537-global-four-arm-emptiness.md)
first proves that global component four-arm labels are empty.  Head
[`e21f6220`](https://github.com/LightChainr/Matching-One/blob/e21f6220775d51450d73d4258c3edf2dbdf0a910/results/p537-finite-collar/REPORT.md)
then lifts the coarse labels to `B_inf(z,1)`, retains outer
`J_B=J_W=1`, sums all corner words, and applies the complete root/Schur plus
axis/tilted P4 projection.  The first and only declared minor—rows
`0->1/1->2`, columns axial2 absent/present—is exactly positive:
`det=+2.6904188461441777e-14`.  This rejects rank one for that declared
coarse aggregation.

That old producer did not close the joint physical collar gate: it
canonicalized the source eight ports separately from the collar, omitted the
simultaneous 12-port incidence/role map, and all selected-cell local `a`, `qa`
and `Ea` sufficient sums were exactly zero.  The determinant therefore could
not yet be named direct local pair-kernel transmission.

The same coarse matrix is strongly stage-nonseparable: endpoint-availability
retention is `60.7191%` for first birth and `2.24641%` for second birth, with
ratio `27.0294`, cross-ratio `chi=0.0369967`, and interaction residual
`+9.9275800e-8`.  These are signed Schur cell-mass ratios, not event
probabilities or a localized source carrier.  The branch-only
[`c958b30b`](https://github.com/LightChainr/Matching-One/blob/c958b30b27abbbb0103f7080688ae3a9b58a24cd/results/p537-finite-collar-joint/REPORT.md)
**paired two-site counterfactual** refinement has evaluated axial2 endpoint off/on on the
same rest background while preserving simultaneous 12-port incidence.  It
exactly reconstructs all four legacy cells and the coarse determinant.  Of
1,390 joint sectors, 118 have an exact nonzero determinant; all 106 pooled
complete rectangles and all 30 rectangles complete in both geometries are
exactly nonzero.  The first full-sector determinant is
`-4.005552609094306e-18`, with an interval width `1e-50`.  Thus the finite
value-level obstruction survives inside common physical sectors and is not a
pure separate-key collision.

The decomposition also changes the interpretation: same-sector determinants
sum to only `5.653307232e-17`, or `0.2101274%` of the coarse determinant;
cross-sector bilinears contribute `99.7898726%`.  The coarse value margins
were already nonzero—row sums
`(+4.3555592e-7,-1.7359553e-7)` and column sums
`(+1.0122295e-7,+1.6073744e-7)`—so an automatic two-defect/extra-arm upgrade
is not licensed.  The exact joint result completes the value gate but does not
itself encode edge adjacency.

The post-closure branch has now supplied that missing decision at
[`22f01e33`](https://github.com/LightChainr/Matching-One/blob/22f01e33c2c386e522b17fc781b3adf70e93548e/results/p537-one-defect-diagonal-edge/REPORT.md).
The frozen full-population-coefficient selected graph contains `135,253` exact Bell/rank transition
classes and `8,250,462` physical pair fibres.  Its stronger kernel-changing
subgraph contains `6,846` classes and `740,950` fibres, with a strictly
nonzero first edge and total signed mass `-4.9488399165e-6`.  The declared
decision is `TWO_INDEPENDENT_DEFECT_GAIN_REJECTED`: first birth contributes
`-5.8210906659e-6` (`117.63%`) and second birth contributes
`+8.7225074943e-7` (`-17.63%` cancellation).  All 12 rank-stage by source-
orbit cells exclude zero.  More sharply, within this fixed-anchor,
nearest-neighbour thermal-`z`, source-present kernel-changing selection, the
no-contact mask contains exactly zero classes and fibres; masks touching either
or both local thermal arms carry all support.  This identifies an exact N25
radius-one contact-local, first-birth-dominated signed four-arm/OPE candidate.
It does not yet give a configuration-level joint-map
witness, annular placement or full original-`U` transmission coefficient.
The `c958b30b` value decomposition and `22f01e33` edge decomposition reuse the
same exact N25 populations and are complementary coordinates, not independent
evidence votes.

The later N25 tensor at
[`df4a64f6`](https://github.com/LightChainr/Matching-One/blob/df4a64f68232eec5aa5b8c8a5d920062aaa7808e/results/p537-one-defect-diagonal-edge/REPORT.md)
keeps contact and birth stage jointly.  With rows `0->1/1->2` and columns
single/double contact it is
`[-2.88380e-6,-2.93729e-6; -5.32257e-6,+6.19482e-6]`, with exact determinant
`-3.3498535471e-11`.  Thus the survivor is not a scalar contact loading times
one scalar birth amplitude.

The frozen held-out N65 production at
[`f46c38c3`](https://github.com/LightChainr/Matching-One/blob/f46c38c3088c1a9f4df8ab0f256b88639f0b34a3/results/p537-contact-stage-n65/REPORT.md)
uses 20M new counter-keyed
samples, four shards, 100 new batches, the independent P45 baseline and no
top-up.  The tensor is
`[-1.57857e-7,-9.22766e-8; -3.09130e-7,+3.69680e-7]`; all four marginal 95%
intervals have the frozen signs, while `Delta_cs` excludes zero at about
`6.14` SE.  The decision is `CONTACT_FUSION_COMPLETION_TRANSMITS`.
Branch-only
[`95a695c7`](https://github.com/LightChainr/Matching-One/blob/95a695c7/results/p537-contact-stage-n65/result.json)
serializes the two independent covariance groups in full and adds a clearly
post-hoc shape test: a single scalar multiple of the N25 tensor gives
`chi2=10.3371/3`, `p=.01591`.  This suggests a running two-coupling shape,
but is not a second prospective vote; both result branches use the same seed
and byte-identical shards.

The same branch then advanced through two dependent reanalyses.  Branch-only
[`bab37f21`](https://github.com/LightChainr/Matching-One/blob/bab37f21/notes/p537-contact-completion-commutator.md)
rewrites the negative determinant as a finite noncommuting completion
commutator and records a zero-parameter two-scale fingerprint.  The fit is
excellent (`p=.95893`), but its `N^-3/N^-29/8` pattern was synthesized after
seeing N25/N65 and remains a mechanism candidate, not prospective field
identification.  Current branch-only head
[`f9ba1ff6`](https://github.com/LightChainr/Matching-One/blob/f9ba1ff6/results/p537-full-t-transport/REPORT.md)
uses the identical sufficient statistics to reconstruct the complete
canonical response: `J65=-.00162250989+/-.00018553008`,
`J65/J25=.29396323+/-.03361398`.  The carrier accounts for `5.892%` of full
`T25` and `2.551%` of full `T65`; moving from the prescribed reference `p` to
the pooled root changes only about `.0155%` of `|J65|`.  This retires carrier
closure as the default explanation.  It does not prove the asymptotic
little-o law, and all N65 analyses remain one dependency block.
An even later
[#537 issue-only gauge audit](https://github.com/LightChainr/Matching-One/issues/537#issuecomment-5490511026)
shows why #275 now comes first: the complete `T_t` is invariant under a common
thermal coordinate shift, but a filtered contact-stage allocation need not be.
Its N9 rational controls and N65 sidecar package have not yet been integrated
into Git, so this is a concrete normalization warning rather than a replacement
score or a reason to overwrite the frozen result.

Parallel closed-unmerged Draft
[#545](https://github.com/LightChainr/Matching-One/pull/545) at
[`e1f19e4c`](https://github.com/LightChainr/Matching-One/blob/e1f19e4c4a5438857e7a68d72243414972eeea85/experiments/p337-landing-minor-exact-20260901/REPORT.md)
adds an all-scale local-collar family and an exhaustive single-axis L4
matching-root control (`det=-2.501463041122436e-6`).  It excludes a stronger
pointwise/common-counterterm identity but omits the paired P4 and full fibre;
it remains an independent algebraic regression rather than another evidence
vote.  The later `c958b30b` and `22f01e33` branches supply the joint value and
edge decisions that #545 did not.  Its
subblock `beta=Q(A)/Q(T)` is not the complete global pooled-root `beta/R`.

Before closure, #544 recorded head
[`e8e9c7cf`](https://github.com/LightChainr/Matching-One/blob/e8e9c7cf42b1d8e41da5c93d3780951e1260fdde/results/p537-aggregate-wedge-l5/REPORT.md)
also completed exact axis L4/L5 populations.  `Psi4` is
`-2.5014630411e-6/-4.0685187142e-7`, while `L^4 chi_perp` is
`-0.6623827/-0.6734050`.  This is compatible with a local `L^-4=N^-2`
thermal-orthogonal correction, but two sizes are not an exponent estimate and
the reduced-block `G4` is not full original-U `T_N`.  After the PR closed, its
same-named branch advanced through the frozen L8 and full-T results and now has
branch-only head `f9ba1ff6`; the reduced-L8 artifact remains pinned at
[`22f01e33`](https://github.com/LightChainr/Matching-One/blob/22f01e33c2c386e522b17fc781b3adf70e93548e/results/p537-aggregate-wedge-l8-mc/REPORT.md).
The 100M L6 block first suggested `L^-9/2`; a separately frozen 100M L8 block
then gave `G4=.0006802041±.0000126970`.  Its interval crosses the preregistered
`L^-4`/`L^-9/2` decision boundary, so the literal verdict is
`UNRESOLVED_MODEL_BOUNDARY` and no top-up is allowed.  The opened L8 result
revealed a post-hoc `L^-35/8` pattern: exact L4/L5 nearly fix its amplitude and
predict L8 within `0.40` SE, but L6 is `3.13` SE low.  This is a reduced-block
normalization hypothesis, not an exponent or full original-U result.  The
size sequence stops here.  Existing root derivatives already give the
normalization-free `Xi_L=N^2 G4 sqrt(P_L)`: exact L4/L5 are
`8.00214/7.99820`, L8 is `8.03977+/-0.15006`, and L6 is
`7.72546+/-0.08469`.  The near-eight normalization is useful context, but its
L6 tension and reduced-observer scope forbid promoting it to an exponent or
running more reduced-`G4` sizes.

**Completed finite predecessor: canonical Kreg's joint Q activation reaches original U beyond nearest neighbours.**
The [completed exact result](results/regular-pair-joint-u/REPORT.md),
`open_pr #267` at `f8e30859f05e86ef35d257fc900f97e74f41e21c`, gives
`J2(25)=∂logQ∂epsilon²U=−.0055194314248394015`, split into
NN `−.001751074454402799` and nonNN `−.0037683569704366022`.
All three exact enclosures exclude zero: global additive first-Q closure
and NN-only global transmission are excluded. This does not exclude all
finite-range contact mechanisms or identify an asymptotic field.
One 4.63-second joint-moment calculation used the same exact N25 population,
with no new independent stochastic evidence.
Execution `410015f5505dc2d8ca0e9ac904f656a4adc9fe86`, `branch_only`,
[reproduces all three values](https://github.com/LightChainr/Matching-One/blob/410015f5505dc2d8ca0e9ac904f656a4adc9fe86/notes/regular-pair-joint-transmission-result.md)
with another complete 2^25-per-geometry implementation and supplies the
adjacent-kernel and coordinate-invariance arguments. This is an
implementation cross-check on the same exact populations, not another
statistical vote.

**The separate noncontact spatial decision is also complete.**
Execution `a237968f1d7a82d26b46e83c58179dbba7f1a908`, `branch_only`,
[rejects the frozen C64/r16 zero null](https://github.com/LightChainr/Matching-One/blob/a237968f1d7a82d26b46e83c58179dbba7f1a908/notes/regular-pair-spatial-transmission-result.md):
`C32/r8=3.6591796875e−5`, `C64/r16=6.85546875e−6`, with the latter's
99% interval `[5.20339728e−6,8.50754022e−6]`. Each block has 200k fresh
independent configurations; its 32 pair readouts per configuration are
correlated. This connected colour contraction and the completed original-U
J2 have different geometries, observers and acquisition blocks.

The exact kernel has `s<=1 shared components => g=0`; at s=2 it is the
signed contrast product `d_x(1)d_y(1)`. Every observed nonzero contribution
in these two blocks has s=2, without excluding rare s=3/4 contributions.
This production result is distinct from the older conditional Gram bound.
[PR #509 at `baa5d33b`](https://github.com/LightChainr/Matching-One/blob/baa5d33b2f87b2868aa0cb9d3f6518c93dbf3bff/experiments/p337-regular-spatial-support-20260901/RESULT.md)
also proves that a nonadjacent pair with at most one
shared exterior occupied component has exactly zero first-Q activation,
while `|Cxy|≤(43/16)Pr(vacant endpoints and at least two shared components)`;
an actual two-component witness has `a=1/16`. This is a support selection
rule, not a numerical bound on original U or a decay exponent.

The earlier `branch_only` [completed N25 activation](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/regular-pair-interaction-result.md)
and [score report](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/results/p337-regular-pair-activation/score/REPORT.md)
give `W_Q=∂logQ∂epsilon U=−.04503611397592696` for fixed
`Kreg=K2+K0`, while its direct Q1 epsilon response is exactly zero.
The N25 question is complete; this mixed response is distinct from the
old positive pure-K2 occupation tangent.

Two [new exact results](results/regular-pair-counterterm/REPORT.md),
`open_pr #267` at `21563da4b0cf721a2aa512901f6ffc966ffa8384`, compress
the remaining mechanisms. The [one-site theorem](notes/regular-one-site-q1-thermal-quotient.md)
rules out retaining a nonzero direct U response in the entire
entry-regular, homogeneous, one-original-binary-site class, including
changes to both vacant and occupied tensors: Q1 is only a common
Bernoulli reparameterization. Separately, the [shared-line activated Gram](notes/regular-pair-counterterm-gram.md)
is `3/2+(alpha−1/2)²/2 ≥ 3/2` for the same real counterterm at both
marks. This is a conditional four-line bound, not global U or a universal
field norm.

**Completed route:** keeping canonical Kreg fixed and changing the acquisition
axis from pair distance to the already defined **thermal/pivotal two-channel
carrier** produced the finite carrier result now summarized above. Closed Issue
[#536](https://github.com/LightChainr/Matching-One/issues/536) preserves an
exact N25 midpoint decomposition (currently issue-only, not a Git artifact),
while open [#537](https://github.com/LightChainr/Matching-One/issues/537) now
owns the P1 nonlocal-remainder and near-critical transport theory:

| source support | observable/rank pivot | kernel pivot | root terms | total `J2` |
|---|---:|---:|---:|---:|
| all pairs | `-.00593584195` | `+.00041410428` | `+.00000230624` | `-.00551943142` |
| NN | `-.00085245745` | `-.00089978682` | `+.00000116982` | `-.00175107445` |
| nonNN | `-.00508338450` | `+.00131389110` | `+.00000113643` | `-.00376835697` |

Thus the first task is no longer to discover whether the q/E and kernel
channels exist. The [bounded carrier preflight](results/p337-thermal-pivotal-preflight/REPORT.md)
has now replayed exactly the first 32 existing counters from each frozen
L32/L64 block: 64 configurations and 5,242,880 pair/site callbacks, with no
new seed. All midpoint identities pass. Only 18 callbacks changed the kernel:
11 at L32 and 7 at L64; L64's nonzero primitives are entirely external
shared-component transitions in shells 2/3, and the finite replay contains no
kernel-preserving topological event. These are descriptive callbacks, not
independent samples or a zero-probability result. They support a sparse,
shell-localized, kernel-change-dominant finite interface, while showing that
the negative N25 observable/rank total is not explained by a dense population
of direct local q/E flips in this subset.

The next task is sharper than a generic Delivery-A tail bound. PR #509 at
`ec3941b0` retains the complete response as the root-conditioned Hessian
`partial_u partial_epsilon Yhat=J_N/A_N`, with `u=M` the matching mean; kernel,
readout, root and slope enter one Schur-projected signed covariance.  The
[corrected finite-transfer contract](notes/p537-signed-landing-transport-contract.md)
allocates the full Schur bilinear to each `(pair,z)` fibre, transports
anonymous component labels under C4, forms axis and tilted laws separately,
then takes P4 value and thermal-jet minors.  The N9 value
`1001/536870912` and N16 `-chi/2` are useful negative controls for a pipeline
that silently drops either midpoint channel; they are not the final
source/thermal minor.

The immediate mechanism decision is now **#275's original-observable and
identifiable-prediction contract**, not another contact production.  First
separate local charged insertions, torus propagating sectors, homology-marked
traces, scalar composites and normalized expectations.  Preserve the original
`q/E`, one physical normalizer and the pooled moving-root `U`; an unnormalized
character zero cannot be promoted to a zero after twist-dependent
normalization.  Then make at least two candidate mechanisms give unit-bearing
prediction vectors on at least two correlated raw coordinates from existing
primitive/C3, rho-child `E_top`, P43/P57 and thermal/modulus assets.  Keep all
allowed amplitudes when checking covariance-weighted design rank.  If the
candidates are not identifiable, prove the rank degeneracy and name the one
missing physical rotation or modulus relation.  This P0 uses existing data and
algebra by default; it does not authorize new Monte Carlo.

#537 remains open P1 theory.  Its exact N25 all-`z` ledger, N85 normalized
tensor and a future N130 `F=C+R` partition are conditional support, not current
production.  Resume a larger-size block only when #275 supplies a frozen
observable/normalizer and opposing prospective predictions, or when an
independently motivated block already contains the four cells.
#537 also gives an all-`L` remote rerouting example:
fixed nearby marks, a pivot at distance `Theta(L)`, unchanged rank/shared-count,
yet `g:1/4 -> 1/2`. This forces a genuine two-scale carrier classification.
Only after this falsifier should the project prove the carrier-resolved
signed/absolute dyadic transmission bound; no counters are added. PR #530's
complete original-U influence, PR #531's eta/xi
quotient and PR #532's exactly-two-shared-component factorization are inputs,
not competing new tasks.  The stopped N100/N400 window is not retuned, and
the old conditional `x=17/4` versus `x=21/4` ratios remain historical
predictions.  [Next Targets](docs/NEXT-TARGETS.md) owns the full portfolio;
further epsilon derivatives or an alpha/counterterm search do not follow
automatically.

In a separate strong-source limit, Draft PR
[#533](https://github.com/LightChainr/Matching-One/pull/533) at `5aa929a6`
has now closed the complete specified axis/tilted sign law under the single
gate `beta=L/m^2 -> 0`. After both exterior and hole clouds, every horizontal
reversal, extra carriers, transverse winding, the companion quotient and
rank-sector odds are included,
`U/A_N=-(L^2/Delta)m^(-(2L+1))[I0(2c)^2-I1(2c)^2](1+o(1))<0`.
The earlier `L^2/m^3` alpha crossover and its nonlocal west-step candidate are
proof artefacts, not remaining mechanisms. Fixed-`m` thermodynamics remains
outside this joint limit. This is a valuable parallel strong-source theorem,
not evidence for the critical macro-window pilot or a substitute for the
q/E-weighted pivotal gate.
The bounded old occupation tangent retains its separate
[`W_N=N V_av`, ratio 2 versus 1](notes/local-pair-size-response-predictions.md)
comparison under its stated assumptions; those predictions do not apply
automatically to the new mixed-Q source.

The [weak-Q path sign separation](results/weak-q-path-comparison/REPORT.md)
and regular-endpoint zero remain intact. Separately, execution's `branch_only`
[uniform bound, `85d5e44b`](https://github.com/LightChainr/Matching-One/blob/85d5e44ba8aed471470373f972c670dc7c82bdcf/notes/closed-source-uniform-projection-tail.md)
now proves **N25 `Ustar<0<Udrop` for every real m≥64**. That finite-coupling
half-line is complete; another m point is not pending.

| Mechanism question | Delivered decision | Attention now |
|---|---|---|
| Does canonical joint Q activation reach original U beyond NN contact? | Yes: J2(25)=−.0055194314248394015 and nonNN=−.0037683569704366022 are exact nonzero; the N25 stage×contact tensor is nonfactorizing, the N65 block transmits its signs, and full-T reconstruction gives `J65=-.00162251+/-.00018553`. | The typed carrier is only `2.551%` of full `T65`, versus `5.892%` at N25, so it is not the default complete mechanism. #537 is P1 theory for the nonlocal remainder/near-critical transport. #275 is P0 for the raw observable, normalizer and distinguishable forward predictions; no new MC by default. |
| Does canonical Kreg transmit a noncontact Q response after occupation averaging? | Yes: C64/r16=6.85546875e−6 and its 99% interval excludes zero. Two fresh 200k blocks; 32 within-configuration pairs remain correlated. | Completed spatial null, distinct from J2. The s=2 signed contrast carries all observed nonzero entries, without ruling out rare s=3/4 events or identifying a field. |
| Does fixed canonical Kreg activate original U through Q? | Yes: W_Q=−.04503611397592696 at N25; its direct epsilon response at Q1 is exactly zero. | Completed first mixed response; the completed joint response above uses the same canonical source, not a fitted counterterm. |
| Can an entry-regular homogeneous one-site completion retain old direct V? | No for the whole original-binary-site class, including both vacant and occupied tensors: only a common Bernoulli parameter survives Q1. | This mechanism is excluded; bounded occupation reweighting, specified singular completions and multi-site vertices are outside the theorem. |
| Can a common regular singlet counterterm remove the two-mark Q interaction? | No in K2+c(Q)K0 with unit K2: the shared-line activated Gram is at least 3/2. | Canonical occupation-summed J2 is now nonzero; its negative directional response is not a consequence of the positive conditional Gram bound. |
| Does the old pure-K2 local tangent reach original U? | Yes: V_av(25)=+.0018155512845251097, and local≠full seam is proved. Its unrenormalized tensor has a physical two-insertion pole. | Separate bounded-occupation size ratio 2 versus 1 under the stated assumptions; do not transfer it to W_Q or reopen a regular one-site rescue. |
| Does the stable colour trace reach original U at Q1? | Yes: B1=−.00190483618 and B1_logQ=+.00503649603, both with zero excluded. | Completed distinct seam response, not an unrun local-to-seam identification or a unique local field. |
| Does the Q4 coefficient determine its generic-Q component by dimension alone? | No for the full physical closure family: multiplicities change `5→6` and `−1→0` in two exact examples. N25 packing proves its packet is already stable. | Preserve the larger-torus boundary; the completed N25 Q1 result needs no further enumeration or colour scan. |
| Does the tied source measure the ordinary site-RC Q tangent? | No at N25: tied Q `+.06308268` versus rank-projected site-RC `−.26982803`, with B control `+.33291071`. | Keep B in matched-size weak-Q work. Regular-endpoint Q activation remains excluded in that endpoint family, not in every torus trace. |
| Is the strong-coupling comparison still waiting for a finite-m witness? | No: the entire N25 real-m≥64 half-line has opposite signs. PR #533 also closes the bounded `L/m` axis capillary sign window with a positive two-gas determinant. | The remaining strong-source boundary is noncompact `L/m`, fixed-m oblique control and wider relative-partition estimates—not another bounded-capillary fit. |
| Can a weighted rank jump alone explain the one-hole response? | No: `Xi_reweight=+4.550327123237` offsets jump `−15.306045530801`, giving total `−10.755718407564073`. Source-independent gain also remains excluded by `R=+27.766563581230237`; hard-endpoint closure survives. | Preserve the full normalized defect operator; reweighting contains both rank types, not a fitted extra source or population share. |
| Is Sstar's global response only its explicit unit q term? | No on the fixed N25 pair: `2V_beta_null=+.07291782829951701`, with zero excluded by exact bounds. | Keep the fixed-coefficient exclusion distinct from arbitrary fitted rank sources and within-K/rank residual claims. |
| Does homogeneous Sstar transmit on the exact N50 pair? | Yes: U=1.0615603876876551 and V=+.0543457826695583; the finite V=0 null is excluded on both complete merged 2^50 populations. | Contract complete. Do not infer N100 scaling or a field identity from this finite response. |
| Does the N25 plaquette result establish larger-N F4 transmission? | No: the independent80M block at N65/85/130/170 is `NOT_EXCLUDED`, with all four intervals still inconclusive against±.5. | Preserve `INCONCLUSIVE_STOP_FIXED_BLOCK_WITHOUT_TOP_UP`; no automatic extra sampling. |
| Does #154's fixed lag1 source provide strong entry/completion/global transmission? | Independent165M block rejects B/C; both net intervals are inside the declared weak band. | This source leaves primary H4 attention; the completed secondary does not rescue it. |
| Do #334's frozen closure and residual-amplitude forecasts transfer? | Two distinct new-population experiments eliminate their specified closure/amplitude models. | Preserve their different estimands; no pooled residual, fitted half-law or automatic extra descriptor. |

[Next mechanism question and priority](docs/NEXT-TARGETS.md) ·
[Claim boundaries](docs/STATUS.md) · [Scientific map](docs/RESEARCH-MAP.md) ·
[Repository-only team handoff](docs/TEAM-COORDINATION.md).
Exploration remains open; priority is not a task lock or permission gate.

<details>
<summary>Earlier discoveries and full repository navigation</summary>

- **Completed local tangent and two-hole obstruction:** execution's `branch_only` [pure-K2 result, `923f66b9`](https://github.com/LightChainr/Matching-One/blob/923f66b979a6b6132875f783106c041ed3c0c1a9/notes/local-four-port-transmission-result.md) gives `V_av=+.0018155512845251097`, with local≠full seam proved. [The L17 four-path obstruction](notes/local-pair-two-insertion-obstruction.md), `open_pr #267` at `5864de49d19898e505e6aecc316b9cb824712c70`, has connected mixed-log residue `1/[2(1+v_x)(1+v_y)]` and Gram residue `1/2`, excluding all-exterior regularity of unrenormalized K2 without proving a homogeneous partition pole. [Crossing recoupling](notes/local-pair-crossing-sector-resolution.md) activates all four H⊗H irreps; bare thermal overlap zero is singlet/standard cancellation, not an RG selection rule. The recorded symbolic run took about .114 seconds, without MC, cloud work or an old-score rerun. The newer no-go and activated-Gram definitions are pinned at `c4e4a821ab8adfc2d628b4cd65b4ce1d52535b56`; the canonical-activation execution snapshot was `2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb`, `branch_only`.
- **Completed stable Q1 seam trace:** [the result](results/n25-stable-colour-q1/REPORT.md), `open_pr #267` at `5c1f9d3b7971a41d07db3c9fa4ac86529c90c199`, gives `B1=∂epsilon U=−.001904836180602413` and `B1_logQ=∂logQ∂epsilon U=+.005036496028411871`, both with zero excluded. The 6.759-second old-count/saved-root reduction added no enumeration, samples, root search or tests. [N25 packing](notes/n25-stable-colour-completion.md) fixes this stable central completion; B1_logQ is a mixed response, not total ∂logQ U. The earlier execution `branch_only` [Q4 score, `54352b2e`](https://github.com/LightChainr/Matching-One/blob/54352b2eefa651ca482ca84837053c792e82c71e/results/p337-s4-trace-transmission/score/score.json) remains `J22=+5.440121494634842e−6`. A common regular-root branch between Q1 and Q4 is unproved, so no crossing is claimed. [Dimension-only continuation counterexamples](results/colour-specialization-gap/REPORT.md) remain a larger-torus boundary, removed for this N25 packet by packing. Earlier execution snapshot `7681eedd938019d977ede41a7d74ee1b88ffbc50` was `branch_only`.
- **Earlier exact endpoint result:** [N25 plaquette transmission](results/decimation-plaquette-u/score/REPORT.md), `b8d043fc493ab6d7f808d0c074571d2fdd8fb60f`, gives `V_F4=+.19441468646090693` and the forced N50 correction `+.5996568681566026`. The [closed checkerboard source](notes/decimation-closed-source-and-global-u.md) is exactly endpoint-invariant. Its complete-source response and fixed cycle/rank split are now completed above; the odd-area, mixed-Smith child pair does not establish a repeated finite decimation or a continuum law.
- **#334's independent1M intervention is complete:** [the result, `d0a9daf1`](https://github.com/LightChainr/Matching-One/blob/d0a9daf1132779205f119e9b4470f4eea9cb89c1/notes/p334-independent-normal-intervention-result.md), with score `1164ba91`, gives T=(3.08520±.39187)e−8. Its3SE interval lies above the frozen δ=1e−8, eliminating complete two-score conditional-label mean closure; the old3.6565e−8 forecast survives. Separately, [#509's600k experiment, `14b2c98e`](https://github.com/LightChainr/Matching-One/blob/14b2c98ed3a252a2fe79ce5e124d9484b23a264f/experiments/p334-prospective-intervention-20260831/REPORT.md), rejects both frozen residual-projection bands. Different sources and estimands keep these as separate decisions, not pooled votes.
- **#154's independent165M decision is complete:** [the official result, `f4999e29`, Draft #509](https://github.com/LightChainr/Matching-One/blob/f4999e29612da16a3650f24d124fb59137f053d7/experiments/p154-prospective-transmission-20260831/REPORT.md) excludes both entry-dominant B and completion-dominant C at N85/N340. All four simultaneous channel intervals lie within±.30, and both net intervals within±.50; W remains `not_excluded`, not an identified theory. This lag1 source no longer receives primary H4 attention: no additional samples, changed lag/source or rescue template follows this result. The [completed M10/M11 secondary, `612df8ec`, `branch_only`](https://github.com/LightChainr/Matching-One/blob/612df8ec1cbe3be3938ee2e1f6183a1aefc6510b/notes/p154-clock-line-secondary-result.md) leaves both lines `not_excluded` on the same fresh data and identifies neither clock. The [completed decisions and next microscopic-to-U target](notes/independent-decisions-final-20260831.md) are now recorded, with [Decision Experiments](docs/DECISION-EXPERIMENTS.md) retaining the fixed rules. The [earlier handoff](notes/independent-decisions-handoff-20260831.md) is pre-score history. No new production or third default P0 is implied.
- **Parallel exploratory result:** [P398 stationary arrangement](results/p398-stationary-arrangement/REPORT.md), Draft #267 `58db48dd`, finds99.7344% static-score compression by primal/dual sizes, yet32/60 profile classes fail exact transition closure. A two-versus-four detach witness motivates incidence predictions; the small arrangement residual supplies10.17% of one integrated response. This finite-model line is not a third default top-priority experiment.
- **Three-team handoff:** [`docs/TEAM-COORDINATION.md`](docs/TEAM-COORDINATION.md) connects 数学研究执行 / 数学研究总览 / 数学研究俯瞰 and preserves actual delivery/dependency states. Routine coordination stays in the repository, without task locks or periodic messages.
- **Earlier spatial discovery:** [conditional primitive-line source response](results/norm4-source-line/REPORT.md) — all six sizes reject the common E-plus-clock prediction inside rank1. [Fixed-K decomposition](results/norm4-source-line-fixed-k/REPORT.md) also resolves spatial source coupling at fixed occupancy and rank. Exact rotation to each torus frame reveals an almost real positive response, about.089–.110: much of the large laboratory-frame angular contrast follows the known frame rotation. This adds spatial information beyond rank population; it does not identify the original global U's lattice-H4 mechanism. [Million-mark endpoint history](results/norm4-source-endpoint-1m/REPORT.md).
- [`docs/REPOSITORY-CONTEXT.md`](docs/REPOSITORY-CONTEXT.md) — dated context recovery and earlier discussion history; this is not a prerequisite for new analysis.
- [`docs/ISSUE-PR-INDEX.md`](docs/ISSUE-PR-INDEX.md) — complete original 464-item snapshot plus dated increments through #497 and an active-PR refresh; historical counts are not live totals.
- [`docs/NEXT-TARGETS.md`](docs/NEXT-TARGETS.md) — highest-information moves and their exact decision outputs.
- [`docs/STATUS.md`](docs/STATUS.md) — authoritative claim ledger, including the boundary between `main` and unmerged frontier work.
- [`docs/RESEARCH-MAP.md`](docs/RESEARCH-MAP.md) — the state/source/observer/geometry/acquisition atlas.
- [`docs/ROADMAP.md`](docs/ROADMAP.md) — default attention order and parallel research portfolio.
- [`docs/HYPOTHESIS-BOARD.md`](docs/HYPOTHESIS-BOARD.md) — bold mechanisms, sharp falsifiers and the minimum shared decision portfolio.
- [`notes/production-priority-refresh-after-factorial-20260830.md`](notes/production-priority-refresh-after-factorial-20260830.md) — historical production-priority handoff, superseded as a queue by the full-context review and current Next Targets.
- [`notes/production-first-model-elimination-20260830.md`](notes/production-first-model-elimination-20260830.md) — the production-first loop and the E_top, P250 and Euler-clock eliminations.
- [`notes/scientific-frontier-delta-20260830-late.md`](notes/scientific-frontier-delta-20260830-late.md) — latest production/exact delta and the mechanism-changing queue.
- [`notes/scientific-frontier-delta-20260830-post267.md`](notes/scientific-frontier-delta-20260830-post267.md) — post-snapshot P250/P321/P333/P334/P337 mechanism changes.
- [`notes/proof-carrying-operator-passport.md`](notes/proof-carrying-operator-passport.md) — proof-carrying elimination and a non-blocking six-coordinate Operator Passport.
- [`analysis/research_ledger.yaml`](analysis/research_ledger.yaml) — machine-readable nodes, sectors, sources, experiments and dependency groups.
- [`results/evidence-ledger/latest.md`](results/evidence-ledger/latest.md) — primary-only predictive evidence view.

</details>

## The current scientific picture

### Completed production is not the same as identified mechanism

**Completed finite source separation:** execution's `branch_only`
[complete-source report, `ec01768f`](https://github.com/LightChainr/Matching-One/blob/ec01768f520e85f1acfd9d3fde9bcf855477254e/results/p337-closed-source-n25/REPORT.md)
gives `V_Sstar=+.12616536341416915` and `V_Bvac=+.33291070842057197`
on Gaussian(5,0)/(4,3), from exact2^25 configurations per geometry.
The [fixed action, `0d19179f`](https://github.com/LightChainr/Matching-One/blob/0d19179f6c6c36fdbb34b2d93e35a9d5fe10dad3/notes/decimation-closed-cluster-gas-action.md)
is `Sstar=C+F4+Bvac=2*beta_null+q−3K+2N+2`. Since common K tilt is a thermal
coordinate change, [the new coefficient-only calculation](results/decimation-cycle-rank/REPORT.md)
gives `V_q=+.05324753511465212` and `V_beta_null=+.036458914149758506`.
The fixed unit-q alias fails by exact bounds; arbitrary fitted rank sources
and fixed-K/rank conditional-mean mediation remain distinct questions.
This .281-second reuse adds no enumeration, root search or random samples.

**Homogeneous N50 Sstar transmission is now exact and complete.** PR #509
at `ef3b2c68` [gives](https://github.com/LightChainr/Matching-One/blob/ef3b2c68f824e29421747c805ea7a505aca41908/experiments/p337-homogeneous-n50-20260831/RESULT.md)
`U=1.0615603876876551` and `V=+.0543457826695583` from the full merged
2^50-state population of each geometry, in about 49.85 CPU seconds.
The finite `V=0` null is excluded. This contract ends here: it does not
automatically request N100 or establish a scaling law or field identity.

**Larger-N F4 remains unresolved:** the independent80M result
[`25ca3635`, `branch_only`](https://github.com/LightChainr/Matching-One/blob/25ca3635ea64655923c32adee4b62d683579cdcd/results/p337-f4-transmission-20260831/scored/REPORT.md)
uses20M/100 batches at each N65/85/130/170. All four simultaneous intervals
contain zero, and none lies wholly inside±.5: `NOT_EXCLUDED` and
`INCONCLUSIVE_STOP_FIXED_BLOCK_WITHOUT_TOP_UP`. N25 positivity does not
settle these production sizes. The [one-hole gain comparison, now completed
at `f5c4a74a`](https://github.com/LightChainr/Matching-One/blob/f5c4a74a20bad8589c39e1034cfb209462110dbe/results/p337-endpoint-defect/score/REPORT.md),
has `U_st=10.755718407564073` against gain prediction`.5313680267777353`:
`R=27.766563581230237>0` and companion `Xi=−10.755718407564073`.
These are dependent exact decisions from one defect block, not separate votes.
The [mechanistic interpretation, `bc17b81d`](https://github.com/LightChainr/Matching-One/blob/bc17b81d502fb1ca3323f5c20f63c544bb31602d/notes/checkerboard-single-defect-global-u-result.md)
keeps the hard-endpoint identity intact: only an alternating child face can
lose ambient rank (`ell≤1`), with `DeltaSstar=3−2k_null−ell`.
The [completed normalized-insertion split](results/defect-reweight/REPORT.md),
`open_pr #267` at `e1b96895`, gives reweighting `+4.550327123236791` and
weighted jump `−15.306045530800864`; their sum is the imported total Xi.
The exact nonzero reweighting enclosure rejects the fixed jump-only model.
Its covariance can contain both rank-preserving and rank-changing
configurations. Only the alternating one-eighth of each old finite
population was enumerated (2^22 per geometry,2.036 seconds overall), with
no new random samples, root search, cloud job or test campaign.
Execution's [same completed split, `9057325d`, `branch_only`](https://github.com/LightChainr/Matching-One/blob/9057325d86d54a8a909ebf4f0d5b15ae7e40f1a5/notes/checkerboard-defect-reweighting-decision.md)
uses `U_st=−Xi`; it is a common-result source, not an independent vote.

**From observed turnover to opposite fixed-law tails:** execution's `branch_only`
[four-coupling score, `a70eeff0`](https://github.com/LightChainr/Matching-One/blob/a70eeff09f51ce2fa0fea5ae637e9191efbf2e1f/results/p337-closed-source-finite-coupling/score/REPORT.md)
finds `U_t(log2)=−1.370778221631` against `U_t(0)=+.126165363414`;
all four prescribed m=2,4,8,16 derivatives have negative exact enclosures.
At least one local maximum lies in `0<t<log2`; uniqueness is not established.
The later `branch_only` [square-family theorem, `762dbaf4`](https://github.com/LightChainr/Matching-One/blob/762dbaf4c3afd9925f7e39b27220274312db4dc4/notes/closed-source-square-family-leading-law.md)
sharpens the old `U→0` bound: `U_star/A_N ~ −(L²−6L+6) lambda^(2L+1)/Delta`.
The [new fixed-comparison result](results/projection-drop-tail/REPORT.md)
instead gives `U_drop/A_N ~ +(L−2) lambda^(2L−2+2/L)/Delta`.
Here lambda=exp(−t), N=L², the companion has ell1≥L+2, and Delta>0.
Deleting only `m^(−r)` breaks the lowest straight-stripe cancellation;
no fitted rank coefficient is introduced. The [colour-gas identity,
`85fd4923`, `branch_only`](https://github.com/LightChainr/Matching-One/blob/85fd492312b597b3fa102ea913e4bcc7aeae2acf/notes/closed-source-local-colour-gas.md)
bounds the pressure-density difference by2t/N, but that fixed-t bound
does not preserve this topological observer. The [proof and size table](notes/topological-projection-reverses-global-u-tail.md)
give N100/N225 predictions without new measurements. The later `branch_only`
[uniform theorem, `85d5e44b`](https://github.com/LightChainr/Matching-One/blob/85d5e44ba8aed471470373f972c670dc7c82bdcf/notes/closed-source-uniform-projection-tail.md)
proves opposite signs for all real m≥64 at N25, without locating the earliest
crossover. At execution snapshot `2690f665bc8029cb2370d3f1efcef5eb2853705c`,
the `branch_only` [Poisson theorem](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-poisson-double-scaling.md)
also closes original pooled-U suppression, including oblique quotients,
when the systole grows and `N/m²→ζ<∞`. This joint limit is not fixed m.
The fixed-m oblique problem is reduced to an [order-25 twist penalty](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-oblique-twist-comparison.md)
and [restricted-sector odds mismatch](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-pooled-sector-odds-bound.md);
neither is controlled by pressure equality alone.

**Newest independent decision, `open_pr #509`, `f4999e29`:** the fixed lag1
conditional policy has net original-U responses+.043470±.043631 at N85 and
+.060675±.082657 at N340. Their simultaneous95% intervals are
[−.071640,.158581] and[−.157394,.278745], both inside the frozen±.50 band;
the four entry/completion intervals all lie inside±.30. Both dominant
numerical templates fail. This bounds the specified source at these two
finite sizes; it neither proves exact zero nor rejects all temporal
transmission. The `branch_only` secondary `612df8ec` has four line-residual
intervals containing zero, so M10/M11 remain unselected and cannot override
the primary decision. The high clock-map gain is sensitivity, not evidence
that this source occupies its high-gain direction.

**Separate P334 independent decisions:** execution's `branch_only`1M source-normal
intervention rejects complete two-score conditional-label mean closure at the
frozen resolution: T±3SE=[1.909579,4.260822]e−8 exceeds δ=1e−8.
This does not reject a first-Jacobian relation alone. The separate `open_pr #509`
600k first-order H-source experiment gives residual ratios.4988857
[.4360616,.5617098] and.5169035 [.4506760,.5831311], rejecting both C0±.25
and C1[.75,1.25]. Those intervals condition on frozen old training point values;
they do not establish a true-population halving or authorize a fitted half-amplitude
replacement. The [exact crosswalk, `ccabada3`](https://github.com/LightChainr/Matching-One/blob/ccabada318b1eeb12ae28d53391b13ab44c116d2/notes/p334-independent-interventions-crosswalk.md)
separates first-score projection coefficients from the source-orthogonal
conditional-mean component; its identities do not make the experiments
independent votes on one field. The [final decision record](notes/independent-decisions-final-20260831.md)
keeps their source definitions, denominators and dependencies explicit.

**Earlier model-space reduction, `open_pr #267`, `3847a5cf`:** the
[fixed-clock comparison](results/p154-fixed-clock-models/REPORT.md) calibrates
only pooled rank1/root responses, then predicts the unchanged directional and
global readouts. The pure-fourfold relative clock fails; scalar compatibility
does not establish a global mechanism. Full paired prediction uncertainty and
the [prospective cost](results/p154-clock-transmission-budget/REPORT.md) are
saved. No new sampling or test suite was run for these two subsecond analyses.

**Earlier norm4 temporal discovery:** [lagged spatial source](results/norm4-lagged-source/REPORT.md), `open_pr #267`, result `4daae57eef5c945aa050a95cd3d5d5d77582161b`. Center `CB+CW` within the previous occupancy/rank, then add one vacant site: the source couples to the next topological hazard despite its imposed same-time global zero. Moving-root rank1 response is−.008817 to−.002971 across six N; the N260/N340 original-U H4 derivatives remain+.37897±.72219 /+1.27653±1.18227. Signed entry/exit kernels are measured, including direct0→2 events. This115.28-second replay marked2.4M old permutations and retained96-dimensional covariance across the original three groups. It is a one-step-before-readout source for each K, not equilibrium final-configuration fugacity or one fixed early injection for every later time.

These old2.4M permutations and P334's old20-batch prefix archive are **discovery data**. Correct paired covariance preserves uncertainty; it does not turn reuse into independent confirmation. All three independent primary decisions above and the declared P154 secondary are now delivered. Current attention goes to preserving their exact exclusions and relating the distinct finite source/observer models; it does not default to more sampling or new features. P154's primary six-coordinate family and secondary four-residual family reuse the same fresh block, not independent evidence. A failed comparison changes the specified model's priority, not the status of its research Issue.

**Earlier equilibrium-source result:** [occupancy/rank projection](results/norm4-global-source-projection/REPORT.md), Draft #267 `8799dfe1`, gives root-motion contributions +.02432–.02478 from occupancy mixing and +.00409–.00458 from fixed-K rank selection. Both reduce root-comoving rank1 population. At N260/N340, H4 source derivatives in the rank-selective component are −6.3361±4.3464 / +10.9424±8.4926; occupancy components are +.3638±1.3019 / +.9170±1.9262. Global and chain comparisons remain unresolved. This .86-second,150-dimensional same-source decomposition required no replay or root finder.

**Earlier spatial output, Draft #267:** [the conditional primitive-line experiment](results/norm4-source-line/REPORT.md) observes spatial source response **within rank1**, using `O4=(vx+i vy)^4/(vx²+vy²)^2` in physical lattice coordinates. The common E-plus-clock model, including its moving-root transport, fails at all six sizes: chi²/4 is3843.85,2596.77,1049.13,1378.11,8479.76,4629.15 for N65,85,130,170,260,340. These paired readouts retain a396-dimensional joint covariance. One112.14-second replay marked2.4M old permutations, with zero new Monte Carlo samples; aggregation took1.31 seconds.

The informative geometric detail is the **post-analysis exact frame rotation**: multiplying the residual by `exp(-4i theta_period)` leaves positive real amplitudes.088750–.110265 across the twelve geometries, with imaginary magnitude at most.0001164. Thus a large laboratory-frame P4 contrast mainly accompanies rotation of an intrinsic line response. The conditional O4 response and the original root/slope-normalized global U are different observers; their connection remains a mechanism question, not a field identity.

The [same-mark fixed-K decomposition](results/norm4-source-line-fixed-k/REPORT.md) now locates a strong part of this association **within fixed K and rank1**, with joint chi²=22534.76/24 and a maximum component of56.27 SE. Even a geometry-specific source depending only on `K,q,E` cannot produce that within-sector covariance. This points to spatial configuration/line reweighting beyond occupancy mixing. It took1.39 seconds and no additional replay, and is a dependent decomposition of the same observations.

The preceding [endpoint1M analysis](results/norm4-source-endpoint-1m/REPORT.md), execution `6bd46ad3`, separates **root motion and finite rank population**. A common microscopic paired-cluster fugacity moves all six matching roots by about **+0.029 per log Q**. After following that root, rank-1 probability decreases by **0.09250–0.09782 per log Q**, resolved at73–160 SE. This excludes a pure common affine-K thermal-clock explanation at these finite roots; it does not identify an energy operator.

The original H4 source derivative is now **−5.9723±4.2736 at N260** and **11.8594±9.1981 at N340**, compared with previous SE19.98/25.87. Both full lineages, **65→130→260 and85→170→340**, retain their original root/slope-normalized U definition. The source-rigid q2/Jordan extensions remain unresolved (nominal p=.48573/.58735), as do generator drift and the earlier three-global-readout common q/E-plus-clock comparisons. This does not undo the new conditional-line rejection: it uses a different spatial observer. Only the endpoints have1M source-marked old permutations; the four cyclic sizes remain100k. The extra900k per endpoint reuse old counters, with zero new random samples, GPU or server operations.

The independent [P40 million-sample original-U response](results/p40-source-thermal/REPORT.md) is also complete: its two-parent density-source response is −.01089±.00889 and−.01234±.01239 (joint nominal p=.28767), while the common-fugacity root slopes are.028985±.000119 and.028905±.000123. The preceding [even-given-odd response](results/p40-even-given-odd/REPORT.md), `56a6267d`, remains a distinct completed readout: raw C is positive at110–138 SE, raw directional C has p=.21156, and its auxiliary geometry-adapted p=.04911 is not independent confirmation.

The preceding [q/source Gram analysis](results/p40-absolute-cluster/REPORT.md), `a4cbf02a`, and [N65/N13020k Phase-E source replay](results/p154-absolute-cluster/REPORT.md), `eb7ef8c9`, remain distinct completed inputs. The [two-phase baseline/source calculation](results/norm4-source-two-phase/REPORT.md) is also finished: high-precision unmarked q/E curves gave no stable source-noise reduction; the remaining uncertainty lies overwhelmingly in source marks. First original-U derivatives, full-chain residuals, source splits, two-phase comparison, conditional-line response and fixed-K decomposition are completed.

**The angular bridge is completed**, `open_pr #509` at `fb01c44a`: soft `W±=R(1±Re[exp(-4i theta_period) O4])/2` components add back original U and its source tangent. All six U− central values are negative, but the total source response remains unresolved. A source centered within fixed K/rank has exactly zero global-U derivative, despite opposite soft-component responses. The projection above locates the globally visible occupancy/rank parts. The subsequent lagged-source discovery and independent165M channel test are both complete; the latter removes this lag1 policy from primary H4 attention at its frozen resolution. These soft weights are not discrete winding classes or a causal birth/lifetime split.

**P334 completed source discoveries:** [the current scientific card](https://github.com/LightChainr/Matching-One/blob/bc0a18c207e3b09f49ea6b6af6601471114d654a/notes/p334-prefix-response-projection-scientific-card.md) retains hierarchy and local A rank, and gives78.20–80.36% four-contact loading shares using new64 responses/old8 clocks on original00: about20% remains. The earlier80–99% used broader receiver-R0. Both are signed loadings, not predictive R². At #509 `04743caf`, held-out J=BG lowers A risk38–39% and center risk54–59%; conditional center-variance response challenges deterministic fixed-prefix translation. These are completed old-prefix discoveries, not a verdict on the independently frozen normal intervention.

**P398 static compression and dynamical obstruction**, `open_pr #267`, result `58db48dd7b6257930ee5e704a37e438c0cf2faf1`: [primal/dual size profiles](results/p398-stationary-arrangement/REPORT.md) explain99.7344% of stationary-score variance in60 classes, yet32 classes have unequal outgoing rates. An exact two-versus-four detach example exposes missing incidence information. The remaining.2656% arrangement variance supplies10.17% of the final+− integrated response because measure/generator terms cancel. The next parallel exploratory prediction concerns size-conditioned primal–dual incidence; this finite model, P334 contact loading and norm4 lagged hazards have different sources/state spaces and do not jointly identify a CFT field.

**Latest context correction:** [N900 is now completed](notes/frontier-increment-20260831f.md), open PR #484 at `5f30397c`: 32M shared counters,800 batches, width variance2.33946±.12039. Both frozen conditional predictions remain compatible (p=.13497/.08418), with shared target/anchor covariance; there is no forced winner. Its old running label is superseded by the report and successful receipts. The [d increment](notes/frontier-increment-20260831d.md) preserves all147 clocks/noise, canonical crossings, P398 current deletion and the symmetric-two-lobe obstruction; the [preceding crosswalk](notes/context-reconciliation-20260831c.md) retains the earlier context recovery.

The [e increment](notes/frontier-increment-20260831e.md) narrows two further physical questions: full birth clocks already determine uniform-blockade **mean** responses, so spatial variance/final-site collisions add information; P398's fully observed instantaneous current still misses reversible geometric propagation, so the first current observer is no longer a task.

The [P418 per-sample archive reanalysis](results/p418-normalized-archive/REPORT.md) is also complete: the old large common masked-spectrum rejection disappears after correcting unequal batch exposures. Common-fit solutions are numerically supported; radius5-only fits and their derived sharing penalties are not reliable. This changes the radius-flow interpretation, not the independent P250 state increment or exact CRT results. Compatibility is not a unique spectrum or field identification.

The 2026-08-31 review read **all 464 Issue/PR bodies, 1,354 discussion comments and seven reviews**. It recovered results that old opening bodies and earlier overview updates had left as “next”: real E_top and #370 production analysis, P218 coalescence, P155 local thermal-null, P40 motif covariance, P255 ordered-filtration proxy, P334 trigger-graph structure, F5 source separation, N112 E_top C3 and W5 periodic gluing. [The context crosswalk](docs/REPOSITORY-CONTEXT.md) distinguishes completed positive, completed inconclusive, and genuinely unrun work. It is a dated recovery, not a recurring audit prerequisite.

**The fixed macroscopic-window comparison is now conditional P1 support, not
highest attention.** With the same canonical source, the single-field
hypotheses predict `T_4N/T_N=2^(−5/4)` for `x=17/4` and `2^(−13/4)` for
`x=21/4`, where `T_N=N²J2_macro`; the two J2_macro ratios carry the additional
factor 1/16. This windowed bilocal projection is not a homogeneous
global-coupling second derivative, an exponent fit or a field identification.
J2(25)=−.0055194314248394015 and its nonNN part exclude the declared global
additive and NN-only closures; they do not eliminate all finite-range contact
mechanisms. Fresh noncontact C64 is a distinct completed spatial decision, and
neither its positivity nor the conditional Gram predicts J2's sign. The old
bounded-occupation ratio remains separate and supplies no mixed-Q scaling law.
It becomes actionable only after #275 freezes the observable, normalizer and
opposing forward predictions; no automatic K3, counterterm fit or Q1-to-Q4
crossing follows. P154/P334/F4 decisions stand; fixed-m oblique twist/odds
control remains parallel. [Next Targets](docs/NEXT-TARGETS.md) owns the queue.

Several positive results now sharpen this choice. P337's F5 even rows separate W_line/JS response vectors at N325/N425 (`chi2/df=149.93/4,246.93/4`), although [P439's matching loading](results/p439-direct-plateau-transport/REPORT.md) remains unresolved. P267's square-bond N112 E_top C3 response is measured and not collinear with the primitive-line readout; it is not a square-site field identification. P437 measures fixed-support high-order topology at 14.97 SE; the earlier noisy estimator is not a general impossibility result. These remain explicitly unmerged source results.

Two analyses completed before the context-review pause are preserved, not rerun: [P267 response-ray comparison](results/p267-response-ray/REPORT.md) rejects amplitude-only transport for A/E/C/W (`19.87177/3`, nominal `p=.00018045`) but not A/E alone; [P334 capacity allocation](results/p334-trigger-capacity-allocation/REPORT.md) separates opposing support and side-imbalance contributions on selected saved graphs. Thus neither first amplitude profiling nor first graph replay is still the next task.

P398 has moved beyond its completed width-four positive matrix and anisotropic spectrum. Draft #267 [`8f7a587`](results/p398-physical-two-point/REPORT.md) gives exact width-five positive-separation Hankel rank8; [`4846adf`](results/p398-fixed-readout/REPORT.md) now measures observer-dependent slow-mode error without rerunning the transfer engine. Branch-only [`552c45d`](https://github.com/LightChainr/Matching-One/blob/552c45d7595ebcb0d04555cec03b2a5bfd8da44a/notes/p398-width8-source-spectrum.md) has already completed continuous width-eight propagation: two symmetry-protected rays, each with 93 propagation directions, and explicit size-two/contact-multiplicity emissions. This is not the width-five discrete C5 character. Protected rays, exact state dimension, useful approximation and continuum fields are distinct; neither first fixed-readout analysis nor first width-eight propagation is pending.

Open [PR #484](https://github.com/LightChainr/Matching-One/pull/484) has now acquired **N100, N400 and N900** in separate2M/8M/32M blocks. N400's common density-map remainder is unresolved (3.90086/6, p=.69009), while its signed odd profile broadens. [`fb1a944`](https://github.com/LightChainr/Matching-One/blob/fb1a944e1ef34e9b9dfcf32c59af25f44ce43d9a/notes/p267-rank-clock-width-decomposition.md) finds96.895%±.507pp of that broadening remains before canonicalization. N900's effective N400→900 width exponent.30688±.03861 is finite-scale, not a new critical exponent. The next model must predict additional area/center/shape or physical-source responses; the first N900 width comparison is complete. Derived views within each block share their covariance.

Two further source updates matter for the same reason. P398 [`c9dc218`](https://github.com/LightChainr/Matching-One/blob/c9dc218f5522502cce8cca539b876ed5faa49b8a/notes/p398-width8-memory-motifs.md) has already decomposed memory into directed motifs: its strong normalized plus response reflects a nearly dark source, not stronger bare feedback; triplet embedding is measured. P334 [`6358ba4`](https://github.com/LightChainr/Matching-One/blob/6358ba49ef390c10a3f501b589ba7ba1d4e05b09/notes/p334-full-physical-birth-clock.md) solves both selected N425 full physical waiting-time laws, and `7401c93` also gives their conditional canonical K2 curves and removable suffix noise. These are completed finite mechanisms, not population H4 attribution. The [context crosswalk](notes/context-reconciliation-20260831c.md) connects them to #398 and #334/#429/#487 rather than assigning their first calculation again.

The Phase-E boundary is now narrower. The integrated score is the derivative alias `P4[S']=P4[E_top']/2`; the measured radius-1 mixed row improves A/E/C by only `1.981<4`. Draft #267 [`c0880c2`](results/p154-fixed-k-interaction/REPORT.md) also completes the local-edge Q/R/H replay on those same 20k samples per size. The fixed-K arrangement response J_R is unresolved (z=.152/.891; joint .8173/2, p=.66455), not evidence that everything is H(K). R contains a global K counterterm, so it is not a strictly local canonical field. The first interaction response is done; the original physical norm-4 residual remains unidentified. The [source-versus-connectivity note](notes/p154-local-source-versus-connectivity-emission.md) explains why this local-source result and P398's connectivity emissions answer different questions.

### One exact finite object, two activations

On honest periodic square-cell tori, the integrated digital Alexander theorem gives, configuration by configuration,

```text
r_black + r_white = 2,
q = r_black - 1 = 1 - r_white.
```

Along a uniformly ordered site permutation, the ambient black homology rank is nondecreasing and activates at most twice. With `K1=K_minus` and `K2=K_plus`,

```text
q_n = -1 + 1{n>=K1} + 1{n>=K2},
M_N(p) = -1 + E[H_K1(p)] + E[H_K2(p)].
```

Thus the matching curve is exactly an equal mixture of first-wrap and second-direction-completion distributions, not generically one latent threshold. The raw midpoint `C_raw=(K1+K2)/2` and raw gap `G_raw=K2-K1` retain information discarded by the scalar mixture; the normalized clock/lifetime coordinates used below are `C=C_raw/(N+1)` and `W=G_raw/(N+1)`.

The main line contains three complementary finite checks. PR #282 identifies `K_minus/K_plus` with the first and second essential-H1 births; PR #283 reconstructs `P0/P1/P2` and the mean rank-one lifetime; PRs #284-#287, #292, #320, #326 and #330 exhaust the short-period quotient frontier through index 13. The cached/subset-DP result covers 140 HNF representatives, 101,140,028,118 weighted filtrations, 285,020 cached subsets and 500,805,335,024 rank-one plateau steps with no birth, reflection, rank-sum, reconstruction or line failures; maximum `iota` remains one with no path evolution. This is exact finite regression and counterexample-search evidence, not by itself an unrestricted theorem.

The `branch_only` parent `73d4960` claims unrestricted rational duality through a canonical four-sheet cover. Its child `c1a72e5` now claims the missing integral saturation: classify honest carrier images on every `qL` cover, then choose a prime coprime to any hypothetical downstairs Smith defect. If that proof survives review, every rank-one component has `iota=1` and the index-13 census becomes a regression alarm rather than inferential evidence for saturation. It remains unmerged and concerns ambient component homology—not graph cyclomatic number, `p_c`, a CFT field or a finite-size law. A separate branch-only theorem proves that the honest-torus balance root `p_L` converges to square-site `p_c`; the open scaling question is now its near-critical window and rate.

The stacked P334 line has now moved beyond a tiny oracle. A real N65 same-modulus 20k-per-orientation archive records `tau1`, the up-to-sign primitive first-born line `ell1`, `tau2` and a separate `DIRECT_RANK2` atom. Exact N13/N17 censuses initially put about `76/24` of the signed response in the two line orbits near `p_ref`, but the later asymmetric-HNF atlas shows that this share and close root spacing are not laws. What survives is sharper: scalar source/sink currents determine their balance times, while the exact character Gram sign determines reinforcement versus cancellation. The `branch_only` chain through `7ef99ae` proves that two-carrier BA plus aggregate TM implies ULC, makes aggregate TM imply every proper Hall cut, and reduces aggregate TM itself to the curvature-corrected Rayleigh inequality `4DF<=M^2+4Y(T-D)`. Later exact configuration gates at `c2d8170 -> 6a7b042` refute universal mark-only, union-preserving and base-only injections: after the fixed-line output-base correction, base-only reaches `588/1152` and transverse-only `768/1152` on the minimal N6 rows. One carrier-base exchange plus one independently released transverse output mark reaches `1152/1152` on every N6 row and `5760/5760` on the first Smith-(2,4) N8 gate. Full Hall is exactly equivalent to a factor-N compressed orbit-Hall problem for arbitrary HNF. The frozen N9 shell sharpens the obstruction: 22/28 rows saturate, while rows 1, 3, 6, 9, 15 and 24 each miss exactly `2160=5M`; their one-mark image is the nonuniform fiber sum `12*0+8*67+4*70+24*164=4752=11M`, and two output marks close target capacity. The live theorem is now a uniform two-mark source injection—or the first arbitrary-HNF coarse min-cut—not another Hall atlas or a return to the refuted switches; BA remains a separate concordance injection.

A `branch_only` exact review at `fee3328` shows the N10 `(k,rank,line)` state is not Markov under uniform-permutation growth: age classes differ in exit hazard by `1/57`. Production `751f8b3` and Draft reuse `2d47d72` localize much clone dependence to checkpoint H2. Saved-graph replay `1b5a9de` and Draft `2e32fd0` then separate scalar collisions, capacity and side organization. Open [PR #491](https://github.com/LightChainr/Matching-One/pull/491), `ab90201`, now explains the sides by cutting an occupied essential cycle: second-rank continuation becomes two-terminal vertex connectivity in the stated embedded rank-one scope. It reconstructs the two N425 witnesses and their W2 difference `540=472+68`. Population directional loading and cut-covariant physical marks remain open; the first cut construction and the selected graph explanation are complete.

P337 gives that archive a finite Fourier basis. Conditional on the unmerged integral-saturation theorem, flat twists by any finite abelian group of order `n` satisfy `S_n=n^2 P0+n P1+P2`; `S2,S3` and normalization invert the entire unmarked rank source. The exact F3 representation is one A4 triplet `[H,A,D]`, not three fields. On Gaussian-ideal quotients the unweighted A/D coordinates are exact quarter-turn-odd nulls; explicit sources `q_A=T_01-T_10` and `q_D=T_12-T_11` activate the missing B1/B2 sectors. The four-generation reuse completes the midpoint split: at all four sizes risk composition opposes the observed completion H4 while conditional hazard has the same sign; N170 resolves the hazard term at `+.001529+/-.000542` (`p=.00607`). The subsequent ambient-line pilot at `4dd62e3` constructs `W_line`, retains more than 98.27% of its norm after JS projection, yet does not resolve the frozen `O_far/O_sep4` source determinant at either N325 or N425. That observer change has since succeeded: branch-only 2d2a9ab uses two F5 D4-even winding-orbit contrasts and resolves distinct W_line/JS response vectors at both sizes (149.93/4 and 246.93/4). The earlier null belongs to its readout, not source equality. Canonical matching loading and cross-geometry transport remain different questions.

The same joint archives also sharpen the collision story. Branch `b887ef3` rejects the correction-free four-size law `D_N=A N^-5/6` (`p=7.46e-4`) while the N340-to-N680 doubling ratio lands on the fixed `5/6` ratio (`p=.951`); `0.8407` is therefore an effective exponent, not a new asymptotic law. Exact branch `1c93bed` proves every direct birth is a one-carrier theta or two-carrier figure-eight event and hence implies six or eight alternating arms only to injectivity scale. Branch `cfb3ead` supplies the reverse bounded local surgery for the correctly **globally typed** events: exterior rank, common-component or two-component partition, relative deck-address determinant and a routable landing word. Its exact local costs are `p^6(1-p)^18` and `p^8(1-p)^16`. Ordinary untyped six arms remain insufficient, and a uniform landing-separation/nondegeneracy probability is still missing, so no arm exponent, amplitude or `N^-5/6` theorem follows. Raw completion winding/lift, basis, transporter, ambiguity and current microscopic geometry remain `not_scoreable`. No historical server was contacted and no N1360 production was started.

Another branch-only bridge identifies one exact finite projector family at Q=1: up to the positive state-sum normalization, the graph polynomial, matching curve and homology balance are all `P2-P0`. The TL open/closed crossing is its controlled infinite-length endpoint, not a finite-rectangle identity. A child oracle makes a sharper conditional prediction: if the homology restriction is Virasoro-transparent and the moving field is ordinary thermal `Q4 epsilon`, then `C_width(rho)/C_width(infinity)=E4(i rho)`, with `C_width=C_N/rho^2`. The frozen fresh N144/576/1296 campaign remains compatible with the fixed scale and conditional E4 shape, but does not identify an amplitude ratio, transparent projector or field. Exact descendants type the finite object as `Tr(TA)-QTr(A)` in a graded crossed/trivial sum and show regular Q differentiation is semisimple. Branch `96df7c8` then makes the source/Gram obstruction explicit at widths 2--4. The minimal standard FK detach lift has a two-modulus affine Hom jet, but branch `ba3135e` proves that endpoint normalization selects `X0=T,V=0` and the nondegenerate widths 3/4 still fail the radical-Gram/source intersection already at zeroth order. A C3 source/landing doublet repairs width 3, but the revealed width-4 C4 charge-one doublet at `7b40ec7` fails exactly at the Gram-self-adjoint gate. Branch `ab3eed81` has now also failed the final one-dimensional C4 charge-two alternating terminal landing character. Scalar, charge-one and charge-two endpoint irreps are exhausted at width 4. The subsequent P398 rooted/landing accumulator and positive propagation already implement that escape. Their memory/motif analysis is also complete; the remaining bridge is physical full-Q realization, scale/geometry response and site-Matching overlap, not another first endpoint mark or rooted interface.

The completed Draft reanalysis scores all ten archived sizes without new Monte Carlo production or a fitted exponent. `K1` is the larger linearized H4 root-shift point estimate at every size, while the magnitude classification is `K1`-dominant at nine sizes and shared at one. `K2` reinforces it at N=65,85,130,145,170,185,290 and has an opposite-sign point estimate at N=265,325,425. The full-curve Bernstein reuse resolves the apparent reversal: the integrated `K2` area has one positive sign at all ten sizes, while each negative root-point value belongs to a local lobe with a delete-one-stable nearby zero branch at approximately `.59462`, `.59635` or `.59694`. The node-minus-`p_bar` offsets are smaller than their position errors, so no significant ordering is claimed. See [`results/two-activation-h4/latest.md`](results/two-activation-h4/latest.md) and [`results/activation-curve-nodes/latest.md`](results/activation-curve-nodes/latest.md); the isolated 16-worker ARM64 replay of the root split is recorded in [`results/two-activation-h4/server-run.json`](results/two-activation-h4/server-run.json).

A further zero-production decomposition writes every activation direction curve as a same-area pooled-density translation plus a strictly zero-area deformation. For K2, the translation at `p_bar` is positive and the deformation negative at all ten sizes; at N265/N325/N425 the deformation alone overturns translation. Thus the local negative values are distribution reshaping, not a reversal of the integrated second-activation direction. The aligned dependency-group covariance is retained in [`results/activation-transport-shape/latest.md`](results/activation-transport-shape/latest.md).

The same complete covariance has now been used for a direct production model-elimination pass. With `A_top=Delta4 F1+Delta4 F2` and `E_top=Delta4 F2-Delta4 F1`, the pure-odd model gives `445.618/10` (`p=1.80e-89`); zero K2- and K1-directional responses give `182.905/10` and `1041.049/10`; the best common projective line gives `28.593/9` (`p=7.59e-4`); and an uncorrected fixed `N^-13/8` E response gives `37.482/9` (`p=2.16e-5`). These are exclusions of the five declared finite parameterizations on one aligned archive—not an exact probability theorem, continuum-field identification or rejection of corrected/multi-field H4. See [`results/etop-production-elimination/latest.md`](results/etop-production-elimination/latest.md).

The P154 norm-4 production archive now has the same canonical audit. Across N65/85/130/170/260/340, both activation contributions have the same sign and all six are K1-dominant; K2 is individually resolved at four sizes and unresolved at N260/N340. `E_top` is negative at every size. The full-covariance distances exclude `E_top=0` (`5324.015/6`), either activation response being zero, one common `E=lambda A` line (`106.665/5`) and one uncorrected `N^-13/8` amplitude (`177.528/5`). Phase D is complete. The first true mixed Phase-E row is also complete and not selected: its `J_bulk` is unresolved at both N65/N130 and improves the A/E/C common-plane score by only `1.981`. See [`results/norm4-two-activation-h4/latest.md`](results/norm4-two-activation-h4/latest.md), [`results/norm4-etop-production-elimination/latest.md`](results/norm4-etop-production-elimination/latest.md) and [`results/p154-phase-e-mixed-plane-pilot/REPORT.md`](results/p154-phase-e-mixed-plane-pilot/REPORT.md).

The global-ray failure has now been converted into a completed two-factor experiment. At fixed determinant 50, the 100k `tau x topology-map` interaction in `(A_top,E_top,C,W)` is `(.0027604,-.0005511,-.00133656,.00088322)` with `chi2=236.756/4`, `p=4.63e-50`; the additive no-interaction model is eliminated. `E_top` alone is only `-.79` SE, and the frozen A+C diagnostic is `2.075` SE, so the result is finite geometry factorization rather than thermal-field identification. The completed Draft a09758e amplitude-only comparison still fails jointly (19.87177/3), although A/E alone is compatible. A third modulus at N50 within Smith-(5,10) is unavailable; the proposed N250 geometry is unacquired and needs its own anchors or a justified scale law. See [`results/p267-etop-tau-topology-factorial/REPORT.md`](results/p267-etop-tau-topology-factorial/REPORT.md).

An additional Draft reuse of the branch-only P205 quotient prism closes the cheapest character question. A six-coordinate GLS selects `H4/H4` for the two activations (`chi2=2.585155/4`, `p=.629455`): both components reinforce, with fitted amplitudes `A1=.508581` and `A2=.295266`; `K2` supplies 36.732% of the signed fitted H4 amplitude. This retrospective small-N result identifies neither a continuum field nor an asymptotic law, but it shows that the ordinary P205 H4 signal is not purely a first-birth effect. See [`results/two-activation-prism/latest.md`](results/two-activation-prism/latest.md).

Each activation contributes a strictly positive beta density for `0<p<1`, so `M_N'(p)>0`; together with the endpoint signs this gives one unique, simple physical root. This is a finite-volume theorem, not a claim that the root has a closed form.

### The global ordinary and charged sectors disagree usefully

Independent square-site primary blocks strongly disfavor global zero and remain compatible with an H4-like ordinary, unmarked response. The completed norm-5 ordinary transfer favors H4 over the tested H12/H8 aliases, although the child block alone does not reject zero.

A different, deck-charged N325 observable instead gives a frozen likelihood ordering of approximately

```text
H8 / H12 / H4 = 71 / 21 / 8.
```

This is not a discovered spin-8 field: H12 and H4 are not both rejected at a strict 5% gate. It is evidence that the charged and ordinary global observables should not be compressed into one scalar angular amplitude.

The completed P250 chain has moved past both support-first and fitted spectra. Draft PR #267 still establishes one useful one-dimensional fact: a common T/A rank-three annihilator closes the old N505 axis window. An independent 80k/400-batch bivariate block preserves signed C4 but rejects common commuting diagonal ranks through three and the frozen Weyl/free-shared-eigenpair rank-five families. The model-free Hankel result at `a770ac9` includes defective/Jordan commuting states: each plus/minus two-charge block rejects rank at most four and is compatible with rank at most five at the frozen `.01` level (`p=.0543/.0655`), while a raw four-channel shared rank at most five fails (`p=2.84e-5`). The old direction-only comparison at `a46ed63` was nonidentifying. An independent preregistered 1.2M radius-five shell at `11130ae` supports both hand-specific extensions, rejects identity-plus-conjugation and leaves only `Alexander_R2_conjugation` above `.01` (`p=.013379`; R3 `.009686`). The new Draft joint-annihilation score at `99d23a7` then restores the discarded residual magnitude and exactly replays the raw shared-rank statistic before scoring maps. It rejects identity and R2 (`p=.005824`) while retaining only Alexander R3 in the primary family (`p=.077486`). The radius-four joint and radius-five direction survivor sets therefore have empty intersection. This is a mechanism-level tension, not two independent votes: the radius-five score reuses the old 80k moments even though its new 1.2M stream is independent.

That next zero-sample object is now complete. One retrospective augmented `40 x 6` operator preserves the old and fresh delete-one influence covariances separately, replays every old residual block and scores all five pre-existing primary maps plus eleven secondary maps. All 16 reject. The result is not bridge-only: the declared plus-hand `20 x 6` old-plus-degree-five rank-at-most-five Schur chart itself rejects (`p=6.01e-7`), so the direction-only hand-extension result does not support a general `5+5` fallback. This closes fixed-map voting on the present blocks but does not identify rank six, a closed transfer algebra, a microscopic quotient or a field. The next discriminator is larger state dimension versus ordered/path-enriched memory. The archive records historical non-designated `Huawei-CodeBuddy-XPk2PZ` provenance and is used here only as committed Git evidence.

The branch-only radius-six Level-S certificate at `33c557b` remains a valid finite-window elimination: endpoint-Hankel `rank<=5,6,7` rejects and `rank<=8` is first compatible; `93eaab1` rejects the full rank-eight R2-plus-conjugation plane. But PR #416 and production result `dbeb29c` give the correct endpoint interpretation: complete spatial autocorrelation has broad positive Fourier support, and the archived window is compatible with that cone. Rank eight is therefore neither a hidden-state count nor the next object to escalate. Adaptive support changes the observer by retaining intermediate interventions and order. Production `0061e4e` plus Draft reanalysis `eb29446` now show that this retained pre-outcome state adds about 31% held-out predictive gain beyond a low-dimensional residue-Fourier baseline. The live unresolved coordinate is the shared `k=36` residual and the absent complete terminal spectrum, not another endpoint rank.

Open PR #385 adds two exact constraints that remain useful after that elimination. An unrestricted ordinary two-mode family can approach a Jordan semigroup inside any finite observation window, so a positive-radius certificate needs an independently justified spectral gap, conditioning/positivity restriction or exact symmetry obstruction. Its conditional length-five C4 theorem gives one such obstruction only after a closed flat Laurent quotient, cyclic source and internal rotation law are established. The theorem survives; the PR's proposed continuation in the formerly selected R2 gauge does not, and `33c557b` supplies no statistical support for a five-state P250 quotient. Apply a character test candidate-independently only if the eventual flat rank and symmetry premises match its theorem.

The ordinary H4/H12 opposite-alias campaign is now complete on `branch_only` head `86e77db`: 600M paired replicas at each of N305/N325 and 300 batches per design. H4-only remains compatible (`2.224/2`, `p=.329`), but zero also remains compatible (`1.711/2`, `p=.425`); `A12=.1474+/-.1911` is unresolved (`p=.441`). This is a clean H12 null, not a new H4 identification. The third alias remains available but has low information under the current matching-odd amplitude; a stronger independently motivated source should precede another expensive alias row.

### “Spin 4” is not a field identity

The main continuum candidates remain

```text
V_(2,+/-2): x=17/4, |spin|=4, four-leg [2] sector
Q4 epsilon:  x=21/4, |spin|=4, thermal singlet descendant/Jordan candidate
```

An unmerged exact global-endpoint calculation finds an ordinary linear selection zero for the four-leg `[2]` carrier under the explicit hypothesis that the matching observable is a regular unlabeled zero-leg generic-Q endpoint. A separate unmerged width-four Q=4 transfer calculation reproduces the singlet-to-`[2,2]` matrix-element zero while showing that charged/twisted traces can still see that block. This narrows the global selection question; it does not prove absence from every trace, singular Q derivative, defect sector or marked observable.

Main-integrated PR #291 adds an exact necessary observable-design gate: on any single C4 direction orbit the nominal spin-4 character is constant and aliases a scalar (`+1` on axes, `-1` on diagonals). Under a declared shared two-component radial/type normalization, axis and diagonal averages give a determinant-`-2` half-sum/half-difference separation; raw orbit-specific amplitudes require calibration, and the gate still aliases higher `4 mod 8` harmonics such as H12. A typed/internal edge law is the other escape. Therefore an “orientation-resolved H4 kernel” must carry multiple calibrated phases or declared internal transformation data; the six seam numerators test colour character, not spatial H4 identity by themselves.

The unmerged double-projector staircase sharpens the ordinary rank-plane question: the vacuum/KdV spin-4 direction lies on the Alexander-even line and has `A_top=0`, while the regular unlabeled `[2]` four-leg carrier also has zero ordinary overlap. Thermal `Q4 epsilon`/Jordan is the first listed candidate allowed by both selectors. “Allowed” is still not “coupled.” A later nine-geometry branch-only production and Cartan decomposition reject the frozen Q4/Q4-Jordan/H4 maps and classify the global line as contact/projective. The completed mean-line N65/130/260 q2 chain and a two-radius local-UV q2 selector also reject their own frozen parameterizations; these are targeted closures, not a rejection of global H4.

The first positive escape is an external bulk Euler observer `O_ext=V-E+F0`. Its far-from-root `J_D4` coupling is enormous at N325/N425 and survives both the radius-2 contact subtraction and direct source-plane projection. The frozen two-observer continuation is complete: N325/N425 each have 2M paired replicas and 100 aligned batches, and the joint determinant score is `chi2_4=4.27897`, `p=.36957`. The current `[O_far,O_sep4] x [J_Dperp,J_S]` basis therefore does not reject one common projective lane.

The same retained batch-by-occupation tables now give an exact source decomposition without new simulation. For pre-insertion occupation count `n`, `E[O_ext|n]=n-2N(n)_2/(N)_2+N(n)_4/(N)_4`; the radius-two nuisance also has an exact conditional mean. This occupancy clock accounts for 53.14% and 56.16% of the raw far-D complex magnitude at N325/N425. After fixed-`n` centering and within-`n` JS projection, a coherent 46.64% and 43.65% remains, with joint `chi2=390067/4`. Thus density-clock and direct-JS explanations are insufficient in this declared metric, but the survivor is not yet a second microscopic field, Q4 identity or exponent. The next source-rank experiment should condition the survivor on endpoint age/completion/collision or add a genuinely new cross-boundary/winding source—not rotate the completed matrix. See [`results/euler-occupancy-clock/latest.md`](results/euler-occupancy-clock/latest.md).

Main now also fixes an operator-mixing identifiability boundary. The four parity-allowed rows have the exact powers `N^-1`, `N^-13/8`, `N^-5/4` and `N^-5/8`, but their structural map has rank four. A two-amplitude prediction requires external dynamics for `f_I1/f_I0` and `f_T1/f_T0`; parity alone supplies neither ratio.

Main also contains an exact Q=1 velocity fingerprint that makes the future Q-tangent test sharper. The four-leg `V_(2,+/-2)` row has `dx/dQ=-5 sqrt(3)/(16 pi)`, while thermal `Q4 epsilon` has `-9 sqrt(3)/(16 pi)`; their velocity gap is `sqrt(3)/(4 pi)`. The generic-loop spin-8 and spin-12 rows in that oracle are controls, not assignments of the experiment-design H8/H12 aliases. The analytic velocities therefore do not need to be rederived. Tiny VJS 1055e221 already closes a finite measure + projector + explicit-field derivative; the missing bridge is its fixed-field scale dependence and overlap with the original site-Matching response.

P333 exposes a further exact ambiguity before that interface can be trusted. The homology lift `W2D-W0D` and critical-polynomial lift `W2D-Q W0D` have the same Q=1 endpoint but normalized first Q-jets differing by `-pi0D`; raw Q tangents therefore depend on the chosen lift/section. Its sibling descriptor records sector weights, normalization, Q-v path, explicit insertion, projector convention and field counterterm. The two heads share one Phase-A L2/L3 bundle, and their counterterm is not yet a Jordan gauge transformation. The next positive control must project the bottom field and compare a scale/modulus log slope under a complete descriptor.

The separate P333 connectivity-radical chain has also crossed its former next target. Join-only generators are commuting idempotents and cannot carry a Jordan block; bare detach/join words create algebraic transients but fail the first-jet Gram condition. At `6c60b0e`, however, `K=D+J-DJ-JD=(D-J)^2` is nonzero rank one, square-zero and H-self-adjoint on the three-mark radical, with an exact isotropic-bottom Jordan chain. The signs make it non-stochastic. The minimal generic-Q detach certificate at `ba3135e` rules out one scalar loop-weight jet; `081a5ed -> e7e6c80` then rules out one mark and the frozen minimal falling-factorial multimark quotients. The first non-scalar repair is exact at width 3, but the width-4 charge-one row fails at `7b40ec7` and the final charge-two endpoint character also fails at `ab3eed81`. Scalar, charge-one and charge-two terminal endpoint irreps are exhausted. The rooted/landing coupled registry is now implemented by P398, followed by positive propagation and memory/motif analysis. Full-Q physical realization and original site-Matching overlap remain separate questions.

The historical endpoint-mark chain is superseded by P398 `5389200`/`afc619c` for formal closure and `e38fe76`/`b35e100` for positive-measure propagation. Width-five `8f7a587` now separates exact rank8 from a useful two-slow-mode description; width-four `dbd4081` fixes the continuous kernel. The current question is their physical interpretation with unchanged AP/landing readouts, not another first matrix. A physical all-Gram/all-Q realization or site-Matching field identity has not thereby been obtained.

### The finite-size state is compact but not scalar

The N145->290 full-curve result rejects a one-multiplier description through a resolved shape direction. Norm-4 subsequently rejects the frozen analytic q2 scalar and common thermal-jet generator. The post-reveal “Jordan plus one conjugation-even mode” recurrence then passes its frozen N520/N680 generation-four pilot at `lambda=1/2` (`1.314/2`, `p=.518` for scalar U; `9.298/10`, `p=.504` for the jet), removing the visible residual tension. Nearly identical `lambda=0` and `lambda=1` scores mean that pure Jordan, the analytic even-mode choice and persistent curvature remain unresolved; the pilot does not establish a nonzero extra mode or a unique transfer matrix. This raises the value of a geometry/operator coordinate chosen to maximize the frozen models' Mahalanobis separation rather than another free lambda fit.

In the annulus line, PR #247 supplies a useful split inside one correlated raw block: the ordinary/local plus-shell control is compatible with one common per-log amplitude (`p=.3564`), the matching-odd minus-shell scalar law fails (`p=.00382`), and matched-cutoff equality is marginally tense (`p=.04544`). A two-state radial recurrence determined from N325/N425 nevertheless passes a held-out N365 third geometry with joint `p=.965`. A newer branch-only sector score makes the boundary sharper: bounded-window rank one is rejected for `A_minus`, named rank-two classes close, and `A_plus` does not require rank two.

The planned same-semantics Gaussian/annulus rectangle is now complete: a shared generator scores `3.93393/6`, while context-specific pairs improve by only `.21622` with bootstrap `p=.40024`. Shared dynamics is not rejected, enrichment has no evidence and the three candidate lambdas remain unresolved. Production coalescence reuse acd38cf also leaves normal/Jordan/generic families underpowered. The later 4daa50c N650 filtration proxy measures an order-sensitive rank-three contrast (43.1466/3), but is not physical Gaussian/annulus AU-UA. The missing physical boundary state and lift are named in the current context map.

### Geometry and acquisition carry different kinds of memory

The exact CRT join says an unmarked final quotient does not remember the order of its factors. An exact path-filtration witness nevertheless shows that intermediate homology ranks can record when rank two first appears. Commit `af7dd01` proves every current endpoint observable factors through word abelianization; `3128e3e` proves fixed-site delete/add overwrites commute or absorb. Open PR #416 proves that the complete uniformly anchor-averaged N505 endpoint autocorrelation has at least 100 nonzero spatial Fourier modes, and branch result `dbeb29c` finds the archived radius-4/5/6 channels compatible with a positive Fourier baseline. The escape is constructive: `6fbbe5e` gives a covariant state-dependent cut/connector, `f54fb8c` finds positive `Rminus` at N325/N425, and `0061e4e` completes N505 at 200k/400. The zero-new-sample Draft reanalysis `eb29446` holds whole batches and shared replicas together across 10 folds. Adding only pre-outcome antisymmetric rank, basis, site-phase and component-change coordinates reduces held-out loss by 30.97%; plus and minus improve 32.87% and 29.08%. This is a direct production-data state increment, not a runner check. It does not compare against the complete terminal spectrum because endpoint fields and periodograms were never stored. The common post-state `k=36` residual suggests a high-frequency landing/orbit coordinate for a future orthogonal readout rather than immediate resampling.

The triangular-site energy/log-pair control has also moved beyond implementation readiness: an open-PR cross-cutoff production score passes its frozen joint model and resolves a nonzero top-partner cutoff shear. Its `kappa_proxy` is gauge-dependent; a top field with two normalized macroscopic spin insertions is the cheapest proposed gauge-invariant continuation.

The branch-only rank-birth algebra supplies a reusable local interface. The exact birth gates `I01` and `I12` give an even rank-increment row `S=I01+I12=M'` and an odd lifetime row `D=I12-I01=-partial_p P1`; every nonzero `D` insertion carries a canonical plateau line `ell`. Complement symmetry separates the normalized midpoint `C` as an odd clock translation from the gap `W` as an even rank-one lifetime. The q-only mean and simple local q2 continuations are closed, the external Euler bridge survives source Gram orthogonalization, and the completed two-observer determinant leaves the current rank-one lane alive. The later F5 readout at 2d2a9ab resolves two source response directions; the earlier far/separated basis was blind to this distinction. Transport to matching and the physical source interpretation remain separate from that positive result.

Two further recovered positive results change this picture. N112 E_top C3 production 2402a33 resolves r1 (142.199/2) and rejects a common child response ray with the primitive-line observer (21.202/2). Fixed-support square-bond production 386db0a and dependent reuse 888af29 measure high-order topology at 14.97 SE, with 99.8485% of that localized energy above degree five. Neither is a square-site energy-field identification. Separately, threshold curves show a reproducible near-collinear tail-deformation direction after rejecting the tested one-term shapes; its amplitude law remains open.

The Gaussian norm-10 commuting square is currently cost-gated, not acquisition-ready. A source-only power branch projects roughly 27.21B and 55.18B samples per N650/N850 lineage for 80% power; 500M would give expected z only about `.47/.34`, and even an ideal fresh covering CRN saves at most 2.34% variance. These are predictions from reused P45/P57 inputs, not target observations. Keep only a small target variance/runtime pilot; do not launch the multi-billion production line.

Main-integrated PR #340 freezes the Issue-27 rational-width prereveal boundary: widths 1-21 and the existing 22-24 forecast are fixed, while the primary 22-24 target artifact is absent. It contains no target score, model ranking or infinite-lattice threshold result. Score once when the primary source arrives, before any refit. This “Stage B” is unrelated to P321's aspect-ratio enrichment stage.

## Five coordinates for every research node

The scientific atlas uses five independent coordinates:

| coordinate | examples | question it prevents us from conflating |
|---|---|---|
| **state** | finite topological state, continuum module, low-dimensional dynamics, microscopic threshold structure | What object is evolving or being identified? |
| **source** | Bernoulli thermal, rank-birth geometric tilt, Potts-Q measure, projector and explicit-field derivatives, relative cluster fugacity, deck charge, boundary source | What perturbation creates the response? |
| **observer** | global unmarked/rank projector, topology-typed, rank-birth line/landing, matching-even/odd, local pivotal, deck-charged, boundary/cross-microscopic | Which sector is actually being measured? |
| **geometry** | Gaussian cover, annulus, modulus/Hecke child, boundary cross-ratio, triangular/square controls | What change supplies the discriminator? |
| **acquisition** | existing-data reuse, exact construction, small pilot, new production | What new cost and independence does the result have? |

The same exponent or harmonic in two different observer sectors is not automatically the same operator. Several views derived from one histogram are coordinates on one random block, not independent evidence.

Semantic adoption is tracked orthogonally to these five scientific coordinates. The historical semantic-inventory snapshot at `main@5ac456d` recorded 59 score paths: 28 direct typed entrypoints, 27 wrapper-covered frozen kernels, 3 migration-required paths and 1 non-applicable utility; the earlier Draft PR #267 snapshot recorded 65 = 28 direct + 27 covered + 9 migration + 1 utility. These are dated inventory snapshots, not current totals after the new scientific scorers. The incorporated main snapshot is `e300609`, including the completed birth/sign/wedge, projection-leakage, sampling-coordinate, terminal-order disk-embedding and partition-gluing controls. Exact support does not add production evidence. Apply the new analyses first; future semantic adapters improve composability without grading evidence or gating work.

PR #339 adds the exact odd-jet ratios `kappa_r=d_r/d_1^r`, their sparse Jacobian and delta-method covariance propagation. PR #349 supplies exact order-five elimination for a synthetic nonlinear coordinate; PR #354 closes a synthetic order-six derivative/parity chain, and PR #351 regenerates the known N9/N8 `(3,2)` Pell polynomials. None is a data estimate or physical assignment. The former real-observable adapter gap is now closed by branch `f5779b9` on eight high-statistics A/E rows. Further method work should be pulled by a new observation or a surviving model class, not placed in front of production.

Merged PR #341 types the existing N130-to-N170 C4-tangent holdout without changing its historical near-zero residual, observer boundary or evidence role. PR #344 does the same for the two Issue43 full-curve kernels; the matching-even and matching-odd summaries remain correlated coordinates of the same frozen runs, not new independent evidence.

PR #342 adds a resource-only gate: the recorded 16-core ARM CPU baseline implies a fivefold end-to-end GPU threshold of about 3.834M paired permutations/s, with exact output equality and deterministic regression required. No qualifying GPU measurement is checked in, so scientific priority is unchanged and CPU remains the present default; this is a measurement boundary, not a task lock.

The next operator assignment is organized as four independent tomography axes: projective birth/completion flux tests **topological boundary mechanism**; a new connectivity/defect tangent tests **source rank beyond the completed rank-one null**; the P321 `E4(i rho)` curve tests **modulus shape**; and calibrated Q4/Z5 seams test **typed charge-spatial selection**. Agreement sharpens an assignment, while disagreement localizes mixing. These axes are not additive votes, especially when they reuse one dependency group or live in different observer sectors.

## Default attention

The single current order is in [Next Targets](docs/NEXT-TARGETS.md). The original norm-4 physical projection remains attention 1; attention 2 now distinguishes effective slow propagation from exact state closure using the completed AP/landing outputs. The other lanes use the already measured factorial, trigger-graph and F5/C3 results instead of repeating their first analyses.

Triangular gauge-free statistics, near-critical homology, Hall theory, the still-unscored PR #230 multi-u response and speculative threshold mechanisms remain visible parallel opportunities. No new simulation is launched by this editorial update, and priorities never become task locks.

## Claim and provenance discipline

The overview uses four lifecycle labels:

- `main_integrated`: present on the current shared `main` line;
- `open_pr`: inspectable in an open PR but not yet integrated;
- `branch_only`: inspectable remote branch, neither promoted to `main` nor treated as canonical;
- `hypothesis`: a proposed mechanism or future discriminator.

Frozen chronology is preserved. Claim-bearing comparisons require identical observable semantics or an exact registered map, and correlated views of one random block are not counted as independent primary evidence. Those are the only evidence constraints; they do not serialize research.

Implementation and reproduction commands are documented in [REPRODUCIBILITY.md](REPRODUCIBILITY.md). This context/priority update does not run or assign another full test suite.

No closed form for square-site `p_c` is claimed. Main now maps why the known exact site-threshold mechanisms do not transfer unchanged: the matching lattice is different, the naive site-to-bond image has overlapping correlations, the natural cell is four-terminal rather than a solved self-dual three-terminal unit, and finite critical-polynomial roots drift with basis. This bounded obstruction does not rule out a new decorated-cell embedding or correlated exact manifold; a candidate must preserve the full law and produce a basis-independent structural consequence. The threshold-origin chain separately fixes block-event semantics, an O(`s^2`) evaluator and a familywise 400-trial/cutoff-373 plan, but contains no fresh IID result. A new methodology-only contract splits threshold distributions into median location, IQR scale and standardized-quantile shape; its committed rows are synthetic fixtures, not collapse evidence. Separately, the microscopic H4 stencil gate excludes only the declared nonnegative axis proxy and points to a negative-phase oblique orbit. Published estimates remain method-specific and are tracked in `data/literature_threshold_sources.json`.

Main-integrated PR #331 is similarly protocol-only: it demonstrates fail-closed exact-rational interval comparison on two hard-coded synthetic artifacts, not a real gadget, certified probability derivation or new bound. PR #384 supplies the exact four-terminal star baseline; PR #386 exhausts the one-internal-vertex structural search; PRs #394--#396 provide all 27 filtered reliability signatures, the bridge/articulation census and bounded planarity (26 planar, K5 nonplanar). Main `b47d99e` now supplies specified-order disk embeddings: 18 of those 27 candidates admit a disk order, and W5 admits only `(0,2,1,3)` up to reversal. Main `eca7d4e` supplies the 15-state join algebra and two-port gluing tensor. W5 remains the only survivor of the separate frozen balance/structure screen, but its scalar root is `.2979305` and branch-only analysis rejects disk self-duality. The spherical-dual terminal projection loses internal-cycle information; the actual disk-relative boundary partition is instead the unique Kreweras complement. Use the existing ordered gluing inputs to construct a concrete W5/face-cycle-with-leaves periodic composition or medial comparison map, not another generic polynomial or gluing tool.

A branch-only exact classification also fixes which polynomial object is worth developing. Under a natural product-local edge-subset contract, the four-terminal site junction is absent from ordinary Tutte/Krushkal-family states, so the matching defect is not a direct named edge-polynomial specialization. It is exactly a derivative of a typed-site rank-image/site-Krushkal quotient with one necessary and sufficient topology source. Main-integrated PR #329 supplies a stable RGS/group-orbit key, and PRs #386/#394--#396 complete the present bounded census, reliability and planarity layer. This redirects exact work toward typed-port recurrence, joint relative-dual boundary connectivity, periodic comparison maps and the connectivity radical; it does not supply `p_c` or a general no-go for correlated gadgets.

## License

MIT. See `LICENSE`.
