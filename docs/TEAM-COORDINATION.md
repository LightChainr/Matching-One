# Three-team scientific handoff / 三队协作入口

Snapshot: **2026-08-31 17:58 CST**. This is a result-routing entry, not a second
priority board or a live process monitor. [Next Targets](NEXT-TARGETS.md) remains
the single attention order; [the machine-readable handoff](../analysis/team_coordination.yaml)
records source pins and the state reported by each team.

## Current complementary work

| Team | Delivered / active contribution | Next handoff |
|---|---|---|
| **数学研究执行** | Completed P334 common-label responses, ensemble Jacobian rank, full thermal curves and their common covariance. Local work; reports no cloud job or tunnel. | E thermal-dipole / age-weighted lifetime moments are proposed, not yet delivered. They would expose redistribution hidden by a small integral. |
| **数学研究总览** | Maintains Draft [PR #267](https://github.com/LightChainr/Matching-One/pull/267), the scientific atlas and this entry. Conditional-line and fixed-K results are complete at `d2a3f445`; now incorporates the execution team's P334 conclusions as cited frontier results. | Connect delivered results to the next discriminator; consume the four packages below without starting duplicate calculations. No server operation. |
| **数学研究俯瞰** | Draft [PR #509](https://github.com/LightChainr/Matching-One/pull/509) is the delivery home for four new complementary calculations, currently **in progress** by team report. | Return question, source/result SHA, report/JSON, aligned batch statistics and compute-release state. An assignment or runner is not a completed measurement. |

The execution result is **branch_only**, on
`analysis/p334-paired-clock-loading-20260831`, with no corresponding PR in the
17:58 query. Its [Issue #334 result card](https://github.com/LightChainr/Matching-One/issues/334#issuecomment-5476674429)
and [PR #267 pointer](https://github.com/LightChainr/Matching-One/pull/267#issuecomment-5476674567)
are discussion links, not code integration. PR #509 was independently checked
open and Draft at `9690d7c1b9bcb19af3fb67b9c0bf5dc16bd5ef12`; its four ongoing
outputs are not claimed present at that head. PR #267 also remains open and Draft.

## Completed P334 science that changes the next question

The [execution synthesis at a5c39de3](https://github.com/LightChainr/Matching-One/blob/a5c39de35bc141f68ca1cc5e68ec134158f4bfec/notes/p334-common-label-response-rank.md)
uses one actual common-label policy, preserving the joint immediate rank / NN-graph
Euler distribution of both geometries. Its source coordinates are
`g±=(L_first±L_second)/2`; readouts are `S=(Y_first+Y_second)/2` and
`D=(Y_first-Y_second)/Δcos4`. The source sign is not an unperturbed H4 field label.

- **Mean and directional source response separate.** At `p_ref=.59274605079`,
  plus→S A is `−5.41618e−5±7.71650e−6` / `−5.87259e−5±1.26600e−5`
  at N325/N425. Minus→D A is `1.53657e−4±3.46721e−5` /
  `1.48394e−4±2.42130e−5`; plus→D remains unresolved.
  In physical orientation coordinates, the mean Jacobian is nearly diagonal,
  with unresolved off-diagonal entries.
- **Ensemble rank is measured; within-prefix rank is a different question.**
  `det(E_prefix J)` for A is `1.22719e−8±4.07190e−9` /
  `1.55912e−8±4.69216e−9`. This supports two mean-response directions and
  challenges a single fixed ensemble response vector. Prefix mixtures can
  produce this rank even when each prefix has one direction. It does not
  determine `E_prefix det(J)` or `E_prefix[(det J)^2]`, let alone a CFT field count.
- **A small E lifetime integral can hide thermal redistribution.** The
  [completed curves](https://github.com/LightChainr/Matching-One/blob/dcd63ace69eebbe80591e332a26830ca70560a85/notes/p334-paired-euler-thermal-response.md)
  give N425 minus→D E `−6.31893e−5±2.09007e−5` near p=.60629 and
  `+4.76565e−5±2.02398e−5` near p=.66784. Main-lobe areas
  `−3.08404e−6` and `+2.60091e−6` nearly cancel; the full integral is
  `−5.25797e−7±2.21120e−6`. The peak locations and pointwise errors are
  exploratory, not selection-adjusted tests; N325 does not resolve the same
  early negative lobe. E here is a topology readout, not an identified energy field.

These are **dependent views of the original 40k paired prefixes total**
(20k at each N, twenty original batches per N), not fresh replications.
The [shared covariance bundle at b582015e](https://github.com/LightChainr/Matching-One/blob/b582015e64e2d8a59e591c4822b14dedaea58b0f/notes/p334-common-label-tangent-joint.md)
already includes 16,948 coordinates per N. First common-policy response,
first ensemble rank, first complete thermal curve and first covariance assembly
are completed inputs. Add genuinely new coordinates with the same batch IDs.

## Four active follow-ups and why they are different

| Package (数学研究俯瞰) | Scientific discriminator | Input / reported compute |
|---|---|---|
| P154 angular bridge | Six-size soft `W±` components of original U and its source derivative, including moving root/slope and exact addback. Connects conditional winding to the original global observer. | Frozen `764595ea`; ZyTrST running/implementing, XPk2PZ spare Ready. Three norm-4 dependency groups. |
| P334 local response rank | Distinct-quartet unbiased det and det² **within the same 00 prefix**. Separates local two-direction response from mixture-generated ensemble rank; det² avoids signed cancellation. | `ffb70969` / `b582015e`; HZsCM6 reported Starting/input preparation. Same original prefixes. |
| P398 linear response | Exact η=0 Fréchet and zero-frequency Poisson response of the finite process. Tests response to a changed generator, beyond the completed ±1/4 finite differences. | `9690d7c1`; TgFr7R running/implementing. Exact calculation, no Monte Carlo evidence count. |
| P334 finite source | Positive common policies at fixed t=±1 with exact joint-class normalization and paired old-tail importance response. Tests nonlinear continuation/saturation beyond an infinitesimal signed histogram. | `0e4db1b8`, `e32a8593`, `959a7fa2`; TV2N0X newly started by that team. No new tails. |

Machine entries are **team-reported allocations at this snapshot**, not our SSH
observations. No other team's jobs, tunnels or keys were inspected or changed.
Current use does not reserve a research question or grant authority to stop a job.
The local workspace index points here; it need not duplicate this scientific queue.

## Lightweight collaboration convention

1. Before a new overlapping source replay or cloud job, send the other teams a
   short question/source-SHA/dependency/compute note. Independent work proceeds;
   this is communication, not approval or a task lock.
2. At a scientific handoff, send **what changed**, exact result commit + report
   and JSON paths, input lineage, original batch IDs/covariance location, and
   actual execution/release state. Mark proposed / implementing / completed
   separately from branch_only / open_pr / main_integrated.
3. Keep one covariance bundle per reused random block. Different observers and
   post-analysis decompositions can discriminate mechanisms without becoming
   independent votes. Return new P334 columns to the existing batch coordinator.
4. 总览 updates the common result map and next question; 执行 and 俯瞰 keep their
   own branches and jobs. No shared-worktree editing, repeated full-repo review,
   automatic merge, Issue lifecycle change or additional permission workflow.

**Near-term synthesis hypothesis:** source parity, prefix mixture and thermal
redistribution may jointly explain apparently conflicting “one-direction” and
“two-direction” summaries. The active local-rank, finite-source and dipole
readouts can separate these possibilities. This is a proposed finite-mechanism
interpretation, not a settled operator assignment.
