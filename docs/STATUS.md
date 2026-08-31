# Project Status and Claim Ledger

**Status date:** 2026-08-31

`docs/STATUS.md` is the **only authoritative current claim ledger**. `docs/ROADMAP.md` ranks information gain; `docs/RESEARCH-MAP.md` maps durable tracks. Branch-local `NEXT-TARGETS`, Issue handoffs and Draft PR bodies are provenance, not competing current status.

## Hard constraints

The repository keeps three claim-bearing constraints:

1. preserve frozen predictions and committed result history;
2. use identical observable semantics, or a registered exact map, for claim-bearing comparisons;
3. do not count correlated views of one raw random block as independent primary evidence.

This cleanup changes navigation and lifecycle only. It does not rewrite a frozen prediction, raw result, historical report, RNG domain or primary evidence block.

## Current decision — completed tests, P0 empty

The #154 165M-permutation test and both distinct #334 independent interventions are complete. The named failed parameterizations are stopped; there is no pending first production, official score or automatic new block. The broader #154/#334/#337 issues remain open at P1; a policy-level stop does not close those mathematical questions.

The exact transmission chain remains `source -> j_in,j_out -> J_grad=j_in-j_out -> original pooled-root/slope-normalized U`. A strong `J_act=j_in+j_out` response alone does not establish global transmission.

The old M0/M1/M2 wording was qualitative planning vocabulary, **not the actual common frozen model family**. #154 tested W/B/C (weak, entry-dominant, completion-dominant); the source-normal #334 test compared complete two-score conditional-mean closure with a fixed positive forecast; the separate four-contact test used C0/C1 residual bands. Those targets and labels are not interchangeable.

**Stop:** no top-up, new lag/source, fifth contact feature or same-block model rescue follows from these results. #337's existing finite theory is the next attention point, not a P0 production assignment. The established C3 finite-size evidence remains valid without unique operator identification.

**Further completed finite calculations:** [homogeneous N50](https://github.com/LightChainr/Matching-One/blob/ef3b2c68f824e29421747c805ea7a505aca41908/experiments/p337-homogeneous-n50-20260831/RESULT.md) now gives original U=1.0615603877 and V_S=+0.0543457827, with full 2^50 counts per geometry and independent derivative review. The finite zero-transmission null is excluded; positive-sign continuation survives without mechanism confirmation. This contract is finished, with no N100/t/epsilon extension. The separate [specified Q1 closed trace](https://github.com/LightChainr/Matching-One/blob/bea717e826df5a22518774b1725ae7bcbe2cb801/notes/p337-q1-closed-trace-transmission-result.md) also has completed finite landing and V=-0.0019048362. The [canonical local pair](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/regular-pair-interaction-result.md) has also completed finite original-U transmission: the direct Q1 response is identically zero, while the fixed mixed response is -0.04503611398. The [new spatial support result](https://github.com/LightChainr/Matching-One/blob/baa5d33b2f87b2868aa0cb9d3f6518c93dbf3bff/experiments/p337-regular-spatial-support-20260901/RESULT.md) excludes first-Q transmission through at most one shared component and bounds the canonical annealed response by (43/16) times the two-component event probability. All 4140 fixed-kernel entries agree between two independent exact algorithms. Remaining work concerns the weighted spatial probabilities/scaling and the same original-U response, not another finite single-mark score. No new production follows.

## Strongest current evidence

| Statement | Level | Current evidence / provenance |
|---|---:|---|
| Square-site matching-odd orientation signal exists | C3 | Independent P43+P57 primary synthesis rejects **global zero**: `chi2=31.1857355515/4`, `p=2.81e-6`; fixed H4 predictions give `3.4622795373/4`, `p=.484`. Canonical archive: `results/server-20260829/P57-norm5-500m`, present on [`5bc84a6`](https://github.com/LightChainr/Matching-One/commit/5bc84a6fa77eea72f41b7c04f66bde683e87165c). |
| Central square-site odd sector is compatible with `DeltaCos4*N^-13/8` | C3 | P31/P32/P37/P43/P50/P57 evidence retained in `analysis/evidence_ledger_manifest.yaml`, present on [`5bc84a6`](https://github.com/LightChainr/Matching-One/commit/5bc84a6fa77eea72f41b7c04f66bde683e87165c). This is compatibility, not unique operator identification. |
| Frozen norm-5 H4 transfer beats H12/H8 aliases | C3 | H4 `0.4163/2`; **H12** `35.1931/2`; H8 `16.0120/2`, same P57 raw block. The **child block alone** does not reject zero (`1.77635/2`). |
| N145->290 full curve is one scalar multiplier | C3 negative | No: three-level transfer `9.3520/2`, `p=.0093`; canonical path `results/server-20260829/P50-n145-n290-fullcurve`, present on [`5bc84a6`](https://github.com/LightChainr/Matching-One/commit/5bc84a6fa77eea72f41b7c04f66bde683e87165c). |
| Pure `P4[S'] ~ N^-5/4` is sufficient | C3 negative | Prospectively falsified; `52.71634/2` on the P48 new geometry, retained in the evidence manifest on [`5bc84a6`](https://github.com/LightChainr/Matching-One/commit/5bc84a6fa77eea72f41b7c04f66bde683e87165c). |
| One **scalar width** explains the higher thermal jet | C2 negative | No: full-covariance norm-5 width score `24.5004/10`; width-corrected q2 `22.2386/10`, retained on [`5bc84a6`](https://github.com/LightChainr/Matching-One/commit/5bc84a6fa77eea72f41b7c04f66bde683e87165c). |
| Rank-2/Jordan is uniquely established | not claimed | No. Frozen norm-4 production rejects the declared scalar q=2 law while Jordan is a borderline survivor; see PR #273 head [`8b26a30`](https://github.com/LightChainr/Matching-One/commit/8b26a30a785bc142a9d17bfed99a8d0e98ddc4dc). Generation 4 leaves lambda 0, 1/2 and 1 nearly indistinguishable; PR #277 head [`3e855ce`](https://github.com/LightChainr/Matching-One/commit/3e855ced4fd98d8979c0b712636b45c2fa54f969). |
| Primitive square-bond spin-4 and square-site thermal-Q4 are the same sector | negative boundary | No. The primitive square-bond response remains a distinct `x≈4` sector; it is not folded into the thermal `x=21/4` candidate. Current archive and exact controls are indexed from [`5bc84a6`](https://github.com/LightChainr/Matching-One/commit/5bc84a6fa77eea72f41b7c04f66bde683e87165c). |

## Exact semantics and controls that remain current

The Issue #43 even-sector channel correction remains

```text
DeltaS_cross = -DeltaS_either
```

with corrected score `0.5700315436/2` and no refit.

Finite **Russo** / chain-rule semantics remain exact:

```text
M'(p) = pivotal_mass_primal(p) + pivotal_mass_matching(1-p).
```

The **pivotal** identity is a semantic control, not a new mechanism vote.

The N=26 frozen finite families remain falsified:

```text
Beta(5,5): first k=5 difference = -96
Beta(7,7): first k=5 difference = +156
```

Small-quotient exact identities and certificates (including N=7/10/13/16/17 controls) may be C5 within their declared finite scope. They do **not** establish a global H4 mechanism, continuum field count or operator identity.

## #334 — two completed independent interventions

The valuable exact core is the projective continuity law

```text
d_p A_ell = j_in,ell - j_out,ell,
J_grad = j_in - j_out,
J_act  = j_in + j_out.
```

This distinguishes population-gradient current from total birth/death activity without claiming either is a continuum field.

Existing discovery is one dependency family, not many independent votes:

- full original paired birth paths: [`9c495ab`](https://github.com/LightChainr/Matching-One/commit/9c495ab13e65f2bc93dc0849ee3b73f88724c4b1), `results/p334-full-birth-archive`;
- exact/selected-checkpoint cut-network theorem: PR #491 head [`ab90201`](https://github.com/LightChainr/Matching-One/commit/ab90201e88409310632812727e0138c56b455644);
- dual-cycle blocker certificates on the same selected N425 examples: PR #492 head [`0e52dba`](https://github.com/LightChainr/Matching-One/commit/0e52dbaeed53dfffa94592e53e38129c179c5078);
- same-stream M/projective-current crosswalk did not establish nonzero M loading: PR #451 head [`bfbceb2`](https://github.com/LightChainr/Matching-One/commit/bfbceb24f4072e5fd2025a2cecb344014adbd9d8). It is not a third live experiment.

The `bc0a18c` freeze has been executed. The [source-normal result, d0a9daf1](https://github.com/LightChainr/Matching-One/blob/d0a9daf1132779205f119e9b4470f4eea9cb89c1/notes/p334-independent-normal-intervention-result.md) uses 1M fresh prefixes: `T=(3.08520 +/- 0.391874)e-8`, with the frozen 3SE lower endpoint above `delta=1e-8`. Complete two-score conditional label-mean closure is rejected; the fixed positive forecast remains compatible. This does not reject a restriction limited to the first Jacobian.

The distinct [four-contact result, 14b2c98e](https://github.com/LightChainr/Matching-One/blob/14b2c98ed3a252a2fe79ce5e124d9484b23a264f/experiments/p334-prospective-intervention-20260831/REPORT.md) uses 300k fresh prefixes per size. The N325/N425 residual ratios are `0.4988857` / `0.5169035`, with intervals `[0.4360616,0.5617098]` / `[0.4506760,0.5831311]`; both exclude C0 and C1. This is one signed-loading projection, conditional on the frozen old predictor/point values, not a full-state sufficiency test. No new half-amplitude model is established.

These are separate new data groups and different targets. They are not two estimates of one “20% residual,” and their related within-block readouts are not additional independent votes. Both blocks are completed branch deliveries; no repeat or extra descriptor is queued.

## #154 — completed lag1-to-U test and policy stop

The established norm-4 result family remains open mechanistically:

- frozen q2/Jordan production: PR #273 head [`8b26a30`](https://github.com/LightChainr/Matching-One/commit/8b26a30a785bc142a9d17bfed99a8d0e98ddc4dc);
- generation-4 target: PR #277 head [`3e855ce`](https://github.com/LightChainr/Matching-One/commit/3e855ced4fd98d8979c0b712636b45c2fa54f969);
- ordinary integrated `J_bulk` is exactly the existing topology-even coordinate, so it is not a second field direction: [`54b3e80`](https://github.com/LightChainr/Matching-One/commit/54b3e80822fa4c407470cd669912c959b9ea4591);
- strong fixed-p source response does not by itself identify global H4 transport: [`56a6267`](https://github.com/LightChainr/Matching-One/commit/56a6267d6a6826a165f93ed3a64a670ca7088180), `results/p40-even-given-odd/REPORT.md`;
- the earlier lagged-source discovery reused old permutations and left its original global `U` response unresolved: [`dd48177`](https://github.com/LightChainr/Matching-One/commit/dd48177340f169c18cd1fc9217101b54090e1e3a), `results/norm4-lagged-source/REPORT.md`.

The subsequent [165M new-permutation result, f4999e29](https://github.com/LightChainr/Matching-One/blob/f4999e29612da16a3650f24d124fb59137f053d7/experiments/p154-prospective-transmission-20260831/REPORT.md) is complete. Both strong B/C templates are excluded at N85/N340; all four component intervals lie inside W's +/-0.30 band. The net intervals `[-0.071640,0.158581]` and `[-0.157394,0.278745]` lie inside the frozen +/-0.50 stop band. W remains `not_excluded`, not an identified physical theory.

Stop prioritizing this particular lag1 policy as the main H4 explanation. This is a bounded weak response, not exact zero or rejection of every temporal source; no sample top-up or replacement lag follows.

## #337 — finite exact frontier, no new production assignment

The [four-profile two-coupling test, f5aa94c3](https://github.com/LightChainr/Matching-One/blob/f5aa94c3c1d619da2272f2409623ed52c876463d/experiments/p337-two-coupling-closure-20260831/results/REPORT.md) has `D3=0.000439154238... > 1/10000`. It rejects common thermal-plus-S profile reparameterization at the fixed N50 root, including p-dependent coefficients. It does not rule out scalar-U-only matching or count continuum fields.

The completed fixed-m precursor results remain at `2690f665`: [uniform rational bounds](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-uniform-projection-tail.md) prove U_star<0<U_drop for the locked N25 pair throughout every real m>=64, with explicit remainders and positive pooled denominators. The [Poisson coexistence theorem](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-poisson-double-scaling.md) separately proves superpolynomial suppression of original pooled-root U for both laws when N/m² tends to finite zeta on growing-systole quotients, including the oblique companion. Neither is a new sampling block or a uniform-in-size fixed-m result.

The unresolved step is fixed-m oblique twist penalty/restricted-sector odds control or an estimator with demonstrated second-moment control at the fixed comparison. The joint-limit contour bound retains an exp[O(N/m²)] bulk factor and does not settle that fixed-m problem; the tiny N25 m>=64 amplitudes do not establish sampling feasibility. Do not reassign the first finite-window proof, fit another source coefficient or restart the stopped P154/P334/F4 blocks; the larger-N F4 fixed block remains inconclusive without top-up. These branch deliveries do not establish a continuum field or homogeneous H4 mechanism.

## Work explicitly stopped as a default

The following loops are not current execution work:

- more replicas of the same N130/N170 local pivotal/tangent rows; **stop adding replicas** to those rows;
- a **third primitive norm-2** generation whose only purpose is another sign flip;
- another **free exponent** fit before a shape/modulus/transmission discriminator;
- rerunning #334 first decomposition, first 147 clocks, first mean-dose or another equivalent current contraction;
- converting a strong local `J_act`, contact or clock response into a claim about global `U` without prospective transport;
- adding new Scientific-card / bold-conjecture handoff blocks to #154/#334 bodies;
- promoting finite exact controls into global mechanism identification;
- treating bounded PSLQ/algebraic exclusions under parked #1 as progress toward a closed form by themselves.

Coalescence, modulus, local pivotal, Q4/Jordan theory, threshold algebra and other exact programs remain useful archived/parked support. They do not constitute an automatic new production queue.

## Explicit non-claims

The project does **not** currently claim:

- a closed form for square-site `p_c`;
- unique global H4 or `13/8` operator identification;
- unique q=2/Jordan/low-rank mechanism;
- that `J_act` drives the original global `U`;
- that tiny finite-quotient C5 controls identify the large-system field;
- that the primitive `x≈4` sector and thermal-Q4 `x=21/4` candidate are the same operator;
- a rigorous new percolation bound.

Historical handoffs and unmerged assets are indexed in `docs/CLEANUP-20260831.md` and `analysis/artifact_registry.yaml`. Claim language in those historical assets does not override this file.
