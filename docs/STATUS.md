# Project Status and Claim Ledger

**Status date:** 2026-08-31

`docs/STATUS.md` is the **only authoritative current claim ledger**. `docs/ROADMAP.md` ranks information gain; `docs/RESEARCH-MAP.md` maps durable tracks. Branch-local `NEXT-TARGETS`, Issue handoffs and Draft PR bodies are provenance, not competing current status.

## Hard constraints

The repository keeps three claim-bearing constraints:

1. preserve frozen predictions and committed result history;
2. use identical observable semantics, or a registered exact map, for claim-bearing comparisons;
3. do not count correlated views of one raw random block as independent primary evidence.

This cleanup changes navigation and lifecycle only. It does not rewrite a frozen prediction, raw result, historical report, RNG domain or primary evidence block.

## Weekly decision — one transmission question

Signal existence is no longer the bottleneck. This week the project attempts **operator/transmission identification once**, rather than opening more coordinates.

Only two issues are live:

- **#334 — independent birth-current intervention.** Test whether a predeclared current/contact model predicts a fresh population-level response.
- **#154 — temporal transmission to the original global observable.** Freeze a quantitative source -> ingress/egress -> `J_grad` -> pooled-root/slope-normalized `U` prediction before any fresh target is read.

The organizing chain is

```text
microscopic source / contact
    -> j_in , j_out
    -> J_grad = j_in - j_out
    -> d_p E_top
    -> original pooled-root / slope-normalized U
```

`J_act = j_in + j_out` is a useful activity diagnostic, but a strong local `J_act` response is **not** evidence that it drives global norm-4 `U`.

The live model family is deliberately nested and small:

- **M0 lifecycle/current-only** — birth/lifetime/current variables are sufficient for the declared target;
- **M1 contact-regulated current** — one frozen contact/source-normal response is required;
- **M2 extra transfer coordinate** — allowed only if defined and frozen **before** the validation block, never as a same-block rescue.

### Hard stop

The next fresh block is a decision/validation block, not another discovery archive. If a frozen model fails, demote or eliminate it; **do not add a feature on the same block to rescue it**. If #154 cannot produce a distinguishable prospective `U` prediction before new data, do not launch the block. No third P0 is opened during this decision.

A failed transmission test is an endpoint, not a request for a thirtieth readout: the current paper-level boundary then remains the established C3 finite-size phenomenology below.

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

## #334 — evidence entering the current decision

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

The next #334 discriminator is already frozen at [`bc0a18c`](https://github.com/LightChainr/Matching-One/commit/bc0a18c207e3b09f49ea6b6af6601471114d654a), `notes/p334-independent-intervention-freeze.md`: a fresh 1M-prefix intervention comparing the declared M0/M1 response with fixed primary `T`, prediction, tolerance and budget. The freeze contains no target result. Historical archives may train and estimate variance but do not enter the independent validation score.

## #154 — evidence entering the current decision

The established norm-4 result family remains open mechanistically:

- frozen q2/Jordan production: PR #273 head [`8b26a30`](https://github.com/LightChainr/Matching-One/commit/8b26a30a785bc142a9d17bfed99a8d0e98ddc4dc);
- generation-4 target: PR #277 head [`3e855ce`](https://github.com/LightChainr/Matching-One/commit/3e855ced4fd98d8979c0b712636b45c2fa54f969);
- ordinary integrated `J_bulk` is exactly the existing topology-even coordinate, so it is not a second field direction: [`54b3e80`](https://github.com/LightChainr/Matching-One/commit/54b3e80822fa4c407470cd669912c959b9ea4591);
- strong fixed-p source response does not by itself identify global H4 transport: [`56a6267`](https://github.com/LightChainr/Matching-One/commit/56a6267d6a6826a165f93ed3a64a670ca7088180), `results/p40-even-given-odd/REPORT.md`;
- latest lagged-source discovery reuses old permutations and resolves temporal entry/exit effects while the original global `U` source response remains unresolved: [`dd48177`](https://github.com/LightChainr/Matching-One/commit/dd48177340f169c18cd1fc9217101b54090e1e3a), `results/norm4-lagged-source/REPORT.md`.

The next #154 output is therefore **not** another source coordinate. Before any fresh block, freeze a quantitative M0/M1/M2 source-to-`J_grad`-to-`U` prediction, equivalence band and joint covariance budget. If the mechanisms cannot make distinguishable `U` predictions, stop rather than sample.

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

Coalescence, modulus, local pivotal, Q4/Jordan theory, threshold algebra and other exact programs remain useful archived/parked support. They are not part of this week’s validation execution.

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
