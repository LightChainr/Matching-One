# Next Targets: Test the Fixed Regular Interaction's Global Transmission

**Updated 2026-09-01.** This is the single attention board, not a permission
system. Parallel work remains open; no Issue is closed or locked.
History and definitions remain in [Decisions](DECISION-EXPERIMENTS.md),
[Status](STATUS.md) and the [scientific ledger](../analysis/research_ledger.yaml).

## Three resolved mechanism questions, one global question left

The canonical regular completion **has been constructed and scored**.
Execution's `branch_only`
[2ba8863f result](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/regular-pair-interaction-result.md)
uses the fixed `Kreg=K2+K0`, with no fitted counterterm. It gives
`∂epsilon U|Q1=0` and the complete original-U activation
**`∂logQ ∂epsilon U|Q1=−0.04503611397592696`** at N25.
This is a different source from the old pure K2 tangent below.

The new [one-site theorem](../notes/regular-one-site-q1-thermal-quotient.md),
`open_pr #267` at `21563da4b0cf721a2aa512901f6ffc966ffa8384`,
closes the entire **entry-regular, homogeneous, one-original-binary-site**
class: even if both occupied and vacant tensors change, the Q1 law only
changes by a common Bernoulli parameter. The moving-root/slope U removes
that parameter exactly. Searching this class for a completion that
preserves the old nonzero direct U response is therefore finished.

The [counterterm-robust two-site result](../notes/regular-pair-counterterm-gram.md)
at the same commit gives

```text
G1 = [[13/8, −1/4], [−1/4, 1/2]],   det(G1)=3/4,
H_alpha = 3/2 + (alpha−1/2)^2/2 >= 3/2,
K_c = K2 + c(Q)K0,   c(Q)=1+alpha(Q−1)+O((Q−1)^2).
```

No uniform real finite counterterm can make this shared-four-line first-Q
interaction additive. Summing the two holes in the contractible L17
exterior divides H_alpha by `(1+v_x)(1+v_y)` and preserves positivity.
This is a conditional physical contraction, **not** a positive global-U
coefficient or a universal field norm. Different counterterms at the two
marks have no such positive lower bound.

| Mechanism statement | Current decision |
|---|---|
| An entry-regular homogeneous one-site tensor retains the old direct Q1 U response | Excluded for the stated whole class, not merely one counterterm. |
| One-site Q activation can be represented by additive independently closed marks in every exterior | Excluded by the uniform positive two-site bound. |
| The irreducible joint activation survives occupation averaging and transmits into global U | Open. The canonical one-insertion mixed response is complete; it does not answer this two-insertion question. |

## Default attention: one fixed global discriminator, one separate size comparison

| Attention | Next mechanism-changing output | Fixed decision and boundary |
|---|---|---|
| Canonical joint transmission | For the unchanged `Kreg=K2+K0`, determine `J2(N)=∂logQ ∂epsilon² U|Q1,epsilon0`, using homogeneous site-average insertion and joint physical contractions. Start with the already defined N25 pair/root; do not scan c'(1). | A model whose first-Q effective log weight is additive and linear in epsilon predicts **J2=0**. Nonzero J2 excludes that global closure. If J2=0, stop claiming that the positive conditional pair interaction necessarily transmits to global U; do not add a fitted counterterm to rescue that claim. The tensor model's unconditional sign is not yet derived. |
| Fixed occupation-tangent scaling | Separately compare **W_N=N V_av(N)** for the old bounded occupation source on `(5k,0)/(4k,3k)` with its k→2k dilation. | [Derived ratio](../notes/local-pair-size-response-predictions.md): R=W_(4N)/W_N=4V_av(4N)/V_av(N). Under the stated single-field/nonzero-loading assumptions, x=17/4 predicts R→2; x=21/4 predicts R→1. These predictions do not automatically apply to Kreg's mixed Q response. N25 is not an established scaling window. |
| Parallel existing work | Named total weak-Q paths with B control; fixed-m oblique geometric twist penalty and restricted-sector odds. | Mixed Q activation is not a total Q-path derivative; the completed finite-m window is not the remaining fixed-m oblique theorem. |

For the primary readout, a site-average tensor at every vacant site has
`epsilon²/N²` times the sum over unordered distinct-site joint closures;
the second derivative includes the factor two. At Q1 all canonical
insertion weights vanish, so a product of separately closed marks would
miss the first-Q joint term. Reuse the original q/E, normalization, pooled
root and slope functional. **Neither Cov(a_x,a_y) nor the number 13/8 is
J2.** A spatial/size analysis follows this fixed estimand, not a new
descriptor search. No production job is started by this attention table;
any new block should carry its stated precision/budget and predictions.
Parallel lines do not require permission from this overview.

## The old finite occupation tangent remains a separate valid object

Execution's `branch_only`
[923f66b9 result](https://github.com/LightChainr/Matching-One/blob/923f66b979a6b6132875f783106c041ed3c0c1a9/notes/local-four-port-transmission-result.md)
gives **V_av(25)=+0.0018155512845251097** for the specified C4 local
four-port tensor. The ports, full q/E covariances, pooled-root movement
and slope terms are delivered. Two-way configuration witnesses prove
that this insertion is **not** the full seam projector. Neither the first
local score nor that identification question remains unrun.

The new [two-insertion calculation](../notes/local-pair-two-insertion-obstruction.md),
`open_pr #267` at `5864de49d19898e505e6aecc316b9cb824712c70`, changes
the next question. In a physical four-path/two-hole exterior,

```text
G(Q)=Tr(Kbar²)=Q(Q−3)(3Q²−9Q+8)/[8(Q−2)(Q−1)],
Res_(Q=1) G=1/2,
Res_(Q=1) partial_epsilon_x partial_epsilon_y log Z_xy
  =1/[2(1+v_x)(1+v_y)].
```

Thus the unrenormalized finite-strength Q1 local tensor family cannot
be regular in every physical exterior. Its already measured **linear**
occupation response remains valid. This is not a claim that the fully
summed homogeneous partition diverges. The [fixed-cut resolution](../notes/local-pair-crossing-sector-resolution.md)
also shows four active colour blocks; bare thermal overlap zero is a
singlet/standard cancellation, not a thermal RG selection theorem.

The bounded occupation reweighting and the joint local tensor contraction
are different nonlinear objects. **Cov(t_x,t_y) is not the tensor's
two-insertion correlator**: both separate marks vanish in the explicit
four-path exterior while the joint tensor closure has the pole. The linear
size comparison uses a different source from the canonical completion.
A numerical ratio decision still needs a stated size window, correction
allowance and fixed precision/budget; no free exponent or fitted mixture
is implied. An order-one W can also reflect modulation of existing
anisotropy, so it is not a unique thermal-Q4 assignment.

## The completed Q4 and stable-Q1 trace results remain inputs

The [Q4 score](https://github.com/LightChainr/Matching-One/blob/54352b2eefa651ca482ca84837053c792e82c71e/results/p337-s4-trace-transmission/score/score.json),
`branch_only`, gives **J22=+5.440121494634842e-6**. Gaussian90° rotation
and the exact rank-normalization factor identify it with our specified
J22; no duplicate enumeration or score is needed.

The new [N25 packing proof](../notes/n25-stable-colour-completion.md)
makes the full stable-Q continuation explicit:
`beta(Q)=I12(Q−3)/(2Q)+I21(Q−3)/2`. I12 means two once-winding essential
components; I21 means one twice-winding component. Existing seam counts
already separate them. The [completed Q1 result](../results/n25-stable-colour-q1/REPORT.md),
`open_pr #267` at `5c1f9d3b7971a41d07db3c9fa4ac86529c90c199`, gives:

| Fixed complete original-U response | Value |
|---|---:|
| B1=∂epsilon U at Q1 | **−0.001904836180602413** |
| ∂logQ B at Q1 | **+0.005036496028411871** |
| B1 from I12 | −0.001945570733316785 |
| B1 from I21 | +0.00004073455271437206 |

Both zero nulls are excluded by rational bounds. This took6.76seconds,
using old counts and the saved root, retaining explicit beta(Q), measure,
pooled-root and slope derivatives. B1 is not the total logQ derivative
of U. A full trace is also not exhausted by the two-essential-cluster term.

The [general specialization counterexamples](../results/colour-specialization-gap/REPORT.md)
remain important for larger tori: full finite-colour projection need not
commute with stable specialization. The N25 bound `c|u|≤2` resolves
that obstacle **for this packet**. Do not reassign the Q1 landing as missing.

The [finite-jet U functional](https://github.com/LightChainr/Matching-One/blob/f43b3674ce29e12629dd790bcbb7370abc5cefbc/notes/closed-source-removable-twist-jet-interface.md)
is already available: a supplied local landing vector must be fully
contracted and thermally differentiated. If using R/(√Q−1), retain
R through quadratic Q order. No extra Q4 seam search, generic certificate
catalogue or fitted colour dimension is required. Q1/Q4 opposite signs
alone do not prove a crossing on a single globally regular root branch.

## Completed controls and windows: do not assign them again

- [Named Q paths](../results/weak-q-path-comparison/REPORT.md): tied-edge
  log-Q response **+.063082681707085**, projected ordinary site-RC
  **−.269828026713487**, prescribed local B difference **+.332910708420572**.
  Same exact N25 population; no new independent evidence.
- [Regular endpoint selection](../notes/weak-q-paths-and-regular-selection.md):
  `ell P_[2](Q)=0` throughout the regular unlabelled one-insertion family,
  including its regular Q derivatives. The completed Q4 and stable-Q1 trace responses
  do not undo this endpoint zero.
- [Uniform N25 sign window](https://github.com/LightChainr/Matching-One/blob/85d5e44ba8aed471470373f972c670dc7c82bdcf/notes/closed-source-uniform-projection-tail.md),
  `branch_only`: **all real m≥64 satisfy U_star<0<U_drop** at each law's
  own original pooled root, with uniform rational remainders. This task is
  complete. Minimum threshold, crossover locations and sampling accessibility
  were not determined.
- [Axis fixed-m winding](https://github.com/LightChainr/Matching-One/blob/575f35ccce850aff3e3120557a2f18475c2e5936/notes/closed-source-fixed-coupling-peierls.md),
  `branch_only`: Sstar, integer m≥256, 4|L, L≥16, uniformly in activity.
  [Poisson joint limit](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-poisson-double-scaling.md)
  also controls the actual oblique pooled U when N/m²→ζ and systole grows.
  It is not a fixed-m oblique theorem.
- The remaining fixed-m pair problem is explicit:
  [twist penalty](https://github.com/LightChainr/Matching-One/blob/70260cc3bb979fc4aa6e0496a71a065f575e9338/notes/closed-source-oblique-twist-comparison.md)
  `limsup Δ_k/k<7τ∞`, plus
  [sector odds](https://github.com/LightChainr/Matching-One/blob/c9d41971f972ceafd7b585c446e1cd93e2d3ebad/notes/closed-source-pooled-sector-odds-bound.md)
  `d≤γL+o(L), γ<2τ` when winding is `e^(−τL+o(L))`.
  Both sources are `branch_only`; pressure agreement alone does not suffice.

## Independent production decisions remain in force

P154's165M fixed lag1 block rejects B/C and meets the net weak band; that
source leaves the main-H4 attention allocation. P334's two fresh interventions
eliminate their specified models with distinct estimands. F4's80M block
remains inconclusive at its fixed budget. None is overturned by the finite
trace theorem or by repeatedly reinterpreting old data. See the
[decision handoff](../notes/independent-decisions-final-20260831.md).
Canonical E_top, first #370 production, K1/K2, defect reweight/jump and the
strong-coupling tails are delivered inputs, not new first-run assignments.

Complex C3, mixed geometry, triangular invariants, boundary Q and
connectivity/defect radical remain parallel exploration. This update absorbed
the completed canonical mixed-U score and added an exact one-site theorem
plus one counterterm-family reduction in0.168seconds. The latter used no
occupation enumeration, random block, new U score or scientific test campaign.
All ten Huawei machines
are authorized through the updated Skill; no cloud action or live availability
check was needed. Coordination stays [in the repo](TEAM-COORDINATION.md).
Draft PR267 remains unmerged. [The older long queue](https://github.com/LightChainr/Matching-One/blob/f405719264c896aa873dd4aae7292795f544ba99/docs/NEXT-TARGETS.md)
is historical, not a second current priority list.
