# Issue / PR increment c: #492–497, #484, #488–490 and #13

**Capture:** 2026-08-31 **14:10:41–14:11:54 Asia/Shanghai** (06:10:41.412–06:11:54.197 UTC).
Main observed: `75fbfec706d5d6dc2210db16a0d41c01ffba6804`.
Local Draft checkout at start: `7cdf1d575f077cebec5cdfb85a25695cc03d1728`.

The [exact JSON snapshot](../analysis/github-context-increment-20260831c.json) preserves original bodies/titles, complete #492–497 and #484 discussions/reviews, full PR head/base/status metadata and all changed paths/blob metadata. These seven PRs have **zero discussion comments, zero formal reviews and zero inline review comments** at capture. #488–490 receive a lifecycle/head refresh only; their earlier full capture is not overwritten. #13's full33-comment body collection is preserved, with its current body and four new/updated comments after the explicitly recorded prior read anchor5473876430 read for this increment. No comment/file-count mismatch was found.

The new-number cutoff is **#497**. Earlier inventories, numbered-item counts, review notes and timestamps remain historical. No new ranked board is introduced.

## Scientific change that affects the next step

### #484: N400 is already acquired and analyzed

Original title remains **“N100 production: shape splitting and zero-area thermal redistribution”**, although its 06:08:02 UTC body update reports new N400 production. Current open/non-draft head is **`ecde7c9132ed35ec1575bd82f11e816722912e6f`**, branch `analysis/etop-modulus-survivors-20260831`.

This is no longer just a N400 design. The changed-path list contains N400 threshold histograms, moments, metadata, acquisition receipts, scores and finite-transport/scale reports. The body reports8M counters per shape and400 aligned batches across three homothetic shapes:24M geometry-pair evaluations but **one8M-counter N400 dependency block**, independent from N100. Periods are all doubled; production was frozen at the earlier `894b3d8` head.

| Completed item | Immutable commit / entry point |
|---|---|
| N400 production | `3e01b495b5b637b0070705e37b4137a9a0ef0d8b`; `results/etop-n400-three-modulus/{REPORT.md,score.json,raw/,logs/}` |
| N100→N400 scale analysis | `b892d9e69bbdb8f5a57dd7d79b951621be4b0522`; [`notes/etop-n400-scale-separation.md`](https://github.com/LightChainr/Matching-One/blob/b892d9e69bbdb8f5a57dd7d79b951621be4b0522/notes/etop-n400-scale-separation.md), `results/etop-n100-n400-scale-transport/{REPORT.md,scale_transport.json}` |
| N400 finite A/E transport | `results/etop-n400-finite-transport-invariants/{REPORT.md,invariants.json}` at the current pinned PR head |
| Rank-clock de-smoothing cited in the updated body | `fb1a944e1ef34e9b9dfcf32c59af25f44ce43d9a`; `results/p267-rank-clock-width/score.json`; body also cites note revision `31255aadbc6d4e03a6fe4cef53eda03f9e2c7f21`. These are referenced results, not automatically main-integrated facts. |

The reported mechanism change is more specific than another secant vote:

- The fixed shared finite A/E density-coordinate discrepancy is now3.901/6, p=.690 at N400, compared with53.914/6 at N100. This is **loss of resolution at current precision**, not proof that the model recovers.
- The fixed-p free-common-secant score1.558/3 and the six-moment cross-scale raw-amplitude profile are already computed. Its amplitude−.0303 has95% interval[−.2068,.1835]; a sign reversal is not resolved.
- The positive observed change is broadening of the full odd clock: its normalized central second moment in the declared z coordinate rises1.440±.0156→2.123±.0738; the signed area fraction within |z|≤2 falls.9527±.000875→.8682±.00615.
- The updated body cites a completed de-smoothing decomposition:96.895%±.507 percentage points of the broadening remains in the rank-clock profile after removing the canonical binomial readout. Its selected two-size effective width exponent.23796±.01352 motivates an N^(-1/4) working fingerprint, **not a newly established critical exponent**.

Do not assign another first N400 run, first N400 secant test, first two-size width comparison or first binomial de-smoothing. Existing outputs now support a conditional third-scale width question and further interpretation of relative E shape; the two-size statistics alone do not establish asymptotics. All derived views of a given size retain that size's covariance/dependency.

### N900: a frozen experiment is published; runtime is not established here

The PR body says N900 is not yet started. Its newer head's subject is “Freeze N900 intrinsic rank-width experiment and acquire two modulus pairs”, and the two added files are a manifest and runner. The immutable [manifest](https://github.com/LightChainr/Matching-One/blob/ecde7c9132ed35ec1575bd82f11e816722912e6f/experiments/etop_n900_rank_width_20260831.json) says `frozen_before_N900_acquisition` and defines:

- Two shapes2i/4i,32M counters,800 batches; the shear shape is intentionally absent.
- Intrinsic signed rank-clock width predictions2.5655353885 versus2.0947508729, including their shared-anchor covariance.
- Source prediction commit `fb1a944e1ef34e9b9dfcf32c59af25f44ce43d9a`.

There are **no N900 raw/result/receipt paths in the captured PR file list**. Thus the capture establishes a published freeze and runner, not completed N900 data and not verified current process state. The title/manifest/body difference is retained rather than resolved by assumption. Do not launch a duplicate N900 acquisition based on the older “not started” sentence; runtime ownership would require a separate bounded check. No runtime or server was inspected here.

### #492: complementary homology carriers on the same two real checkpoints

[#492](https://github.com/LightChainr/Matching-One/pull/492) supplies a dual-cycle blocker certificate at `0e52dbaeed53dfffa94592e53e38129c179c5078`, branch `research/p429-dual-cycle-blocker-20260831`. Its actual inputs are the existing two N425 checkpoints from `6147e22f53902a94e5f133739f2c1d423691d0b8`, later studied at `c827cd8/1b5a9de`; this is not another random block.

It records two vertex-disjoint complementary-white essential cycles per witness:

| Counter | Cycle lengths | Nonisolated side sizes | Pair edges | W2 |
|---|---:|---:|---:|---:|
|43042514269|20,25|12,14|108|926|
|43042505280|19,43|5,29|108|1466|

Cycle windings are±(12,−19) in period matrix[[425,268],[0,1]]. The body reports29,756 safe-pair comparisons, retained ordered vertices/masks, an independent potential-union-find verifier and exact reproduction of the old graph quantities. Those are source-reported checks; this increment did not rerun them.

The finite statement is conditional: if two essential white cycles overlap only at forced singleton-trigger sites, minimal safe pairs cross their disjoint safe portions. Existence of such a cycle pair for every checkpoint/HNF is **not** proved, and search failure is not automatically an odd-cycle witness.

Relative to the already delivered #491 cut network, this is a **complementary carrier certificate**, not a replacement for the full update-closed vertex-connectivity state. It completes the carrier-position explanation for these same two examples, but does **not** finish population H4 loading of incidence/overlap or longer-horizon vertex reliability. Also do not reopen #491's embedded-NN bipartiteness theorem merely because this different sufficient two-cycle condition lacks a general packing proof. Extending the cycle certificates to the already selected22 checkpoints is a concrete exact-geometry option; a packing obstruction would be a result about that carrier condition, not automatic failure of the cut-network representation.

Relevant pinned paths: [proof note](https://github.com/LightChainr/Matching-One/blob/0e52dbaeed53dfffa94592e53e38129c179c5078/notes/p429-dual-cycle-blocker.md), `results/p429-dual-cycle-blocker/{certificate,verification,metadata}.json`, `scripts/verify_p429_dual_cycle_blocker.py`. All eight changed paths are preserved in the JSON.

## Six new PRs: exact pins, original titles and paths

All six are open and non-draft at this bounded capture; their code and scientific maturity are separate from merge decisions.

| PR | Original title | State | Full head commit | Head branch | Entry path |
|---|---|---|---|---|---|
| [#492](https://github.com/LightChainr/Matching-One/pull/492) | P429: certify dual-cycle blockers on the two real N425 trigger graphs | open, non-draft | `0e52dbaeed53dfffa94592e53e38129c179c5078` | `research/p429-dual-cycle-blocker-20260831` | [`notes/p429-dual-cycle-blocker.md`](https://github.com/LightChainr/Matching-One/blob/0e52dbaeed53dfffa94592e53e38129c179c5078/notes/p429-dual-cycle-blocker.md) |
| [#493](https://github.com/LightChainr/Matching-One/pull/493) | Enumerate all local serial partition monoids | open, non-draft | `234a733ce4a045b032fc0a4118ad04cc8337f176` | `exact/issue-13-serial-local-monoids` | [`analysis/terminal_partition_serial_local_monoids_certificate.json`](https://github.com/LightChainr/Matching-One/blob/234a733ce4a045b032fc0a4118ad04cc8337f176/analysis/terminal_partition_serial_local_monoids_certificate.json) |
| [#494](https://github.com/LightChainr/Matching-One/pull/494) | Certify serial partition power dynamics | open, non-draft | `1880579e1bea821051d23b59f8af23014e13b306` | `exact/issue-13-serial-power-dynamics` | [`analysis/terminal_partition_serial_powers_certificate.json`](https://github.com/LightChainr/Matching-One/blob/1880579e1bea821051d23b59f8af23014e13b306/analysis/terminal_partition_serial_powers_certificate.json) |
| [#495](https://github.com/LightChainr/Matching-One/pull/495) | Certify the idempotent-generated serial sector | open, non-draft | `221ee164bdedb2bd43a266b73683d951ef3c4114` | `exact/issue-13-serial-idempotent-sector` | [`analysis/terminal_partition_serial_idempotents_certificate.json`](https://github.com/LightChainr/Matching-One/blob/221ee164bdedb2bd43a266b73683d951ef3c4114/analysis/terminal_partition_serial_idempotents_certificate.json) |
| [#496](https://github.com/LightChainr/Matching-One/pull/496) | Certify generalized inverses of serial partitions | open, non-draft | `19c38451444402dae5319052862b8c95b734e578` | `exact/issue-13-serial-inverse-census` | [`analysis/terminal_partition_serial_inverses_certificate.json`](https://github.com/LightChainr/Matching-One/blob/19c38451444402dae5319052862b8c95b734e578/analysis/terminal_partition_serial_inverses_certificate.json) |
| [#497](https://github.com/LightChainr/Matching-One/pull/497) | Certify serial stabilizers and cancellation | open, non-draft | `ecc815ad5db02a652d1fb45cb8f1dc2ba7660d40` | `exact/issue-13-serial-stabilizer-census` | [`analysis/terminal_partition_serial_stabilizers_certificate.json`](https://github.com/LightChainr/Matching-One/blob/ecc815ad5db02a652d1fb45cb8f1dc2ba7660d40/analysis/terminal_partition_serial_stabilizers_certificate.json) |

#493–497 consume the same frozen15-state typed serial table: local eSe monoids for12 idempotents; power profiles separating12 idempotents, two index-two collapses and one period-two unit; the14-state idempotent-generated sector; generalized inverses/units; and stabilizer/cancellation tables. Their finite results are now available in open PRs. They do not constitute five new reliability/threshold measurements or prerequisites for the physical continuation analysis.

## Lifecycle increment only: #488–490 are now merged

Their old open states in increment b remain historical. This refresh reads the actual merge metadata; no reviews/tests were re-read or run for these three.

| PR | Original title | Original full head | Actual merge commit | Merged at UTC |
|---|---|---|---|---|
| [#488](https://github.com/LightChainr/Matching-One/pull/488) | Exact: enumerate typed serial monoid congruences | `6c0ed6e7c584e8eb619bd4e50d66bf23abe28471` | `3b2bdb851ecc241548561a6454237d38ed4920da` | 2026-08-31T05:45:31Z |
| [#489](https://github.com/LightChainr/Matching-One/pull/489) | Exact: freeze canonical shortest serial words | `3e549e112964b68f8d60e01a81421b5336c59405` | `60b2877418b311577dd0945aa6f22ed09fd3193c` | 2026-08-31T05:46:05Z |
| [#490](https://github.com/LightChainr/Matching-One/pull/490) | Exact: enumerate typed serial monoid endomorphisms | `e15c1e2d94b4a1e93847b66aa80767b34b06e93a` | `75fbfec706d5d6dc2210db16a0d41c01ffba6804` | 2026-08-31T05:48:05Z |

## #13 discussion increment: source material, not execution instructions

Original title: **“[P2] Automated self-dual gadget and critical-manifold search”**. The issue remains open/unlocked and assigned to LightChainr. Its present body preserves the original proposal and says that existing typed-serial construction does not itself supply a probability comparison or threshold formula.

The prior explicit read anchor is [comment5473876430](https://github.com/LightChainr/Matching-One/issues/13#issuecomment-5473876430), recorded in the existing index. Four later comments were read in full:

| Comment | Created at UTC | Current interpretation |
|---|---|---|
| [5474024747](https://github.com/LightChainr/Matching-One/issues/13#issuecomment-5474024747) | 2026-08-31T05:12:15Z | Completion report for #479–483: one frozen-table result family; included numbers are preserved, not recomputed. |
| [5474072070](https://github.com/LightChainr/Matching-One/issues/13#issuecomment-5474072070) | 2026-08-31T05:18:22Z | Claim for congruence/endomorphism/shortlex work, now completed by the merged #488–490; not a current unfulfilled assignment. |
| [5474291194](https://github.com/LightChainr/Matching-One/issues/13#issuecomment-5474291194) | 2026-08-31T05:48:43Z | Reports #488–490 merged results: 13 congruences, 120 shortlex target records, 71 unit-preserving endomorphisms and 38 idempotent retractions. |
| [5474368154](https://github.com/LightChainr/Matching-One/issues/13#issuecomment-5474368154) | 2026-08-31T05:58:53Z | Claims five further finite-table tasks, now delivered in the open #493–497. Its references to previously ongoing work do not override current merge metadata. |

The completion report's wording that periodic gluing remains missing is not a reason to discard the already completed #438 periodic pair. Keep the distinction between construction of that concrete pair and an unproved probability comparison/self-dual critical equation. Task claims, CI requirements and resource requests inside these comments were preserved as authored context, not executed or adopted as new permission rules.

## Scope and handoff

The largest actual next-step change is already-completed N400 production plus width/de-smoothing outputs. The #492 addition enriches the same-checkpoint geometric mechanism; the five new monoid PRs are finite support work. #488–490 require lifecycle refresh, not rediscovery.

Only this report and its JSON snapshot were written. No INDEX/ledger/registry edit, scientific calculation, test, Monte Carlo, GitHub mutation, commit, push, merge or server/process inspection occurred. No item above#497 was fetched.

## Later bounded overlay: lifecycle and active discussions

This is a **separate 14:16:43–14:21:03 Asia/Shanghai read** (06:16:43–06:21:03 UTC), not a correction of the original14:10–14:11 snapshot. PR/Issue reads stopped at06:19:59.663 UTC; the final main-ref-only observation was06:21:03.988 UTC. New-number cutoff remains#497. The JSON `bounded_overlay` preserves exact additional bodies/comments, the all-open PR head list and the non-atomic capture times.

### Actual completion since the first capture

#493–497 were merged after the initial open snapshot. Their full original heads remain in the first table; actual merge metadata is:

| PR | Merge commit | Merged at UTC |
|---|---|---|
| #493 | `621368b39e28a1ec80aad4d80fcf4fc49bc1f638` | 2026-08-31T06:14:16Z |
| #494 | `d777be4508458598ef2853bb66b458a2bb9e9ee5` | 2026-08-31T06:14:19Z |
| #495 | `ba1991ec05cc66fcd06ae3d57666172eef4506ca` | 2026-08-31T06:14:21Z |
| #496 | `331ba6d17ce7acadbbfe97b5648d03b5c0fca513` | 2026-08-31T06:14:23Z |
| #497 | `2ee16678b64a065d6d235ef69941c8a46f3a3878` | 2026-08-31T06:14:26Z |

The observed GitHub main ref is now `2ee16678b64a065d6d235ef69941c8a46f3a3878`. This is a remote metadata observation, **not a fetch or merge into the Draft checkout**. [Issue13 comment5474492110](https://github.com/LightChainr/Matching-One/issues/13#issuecomment-5474492110),06:14:45 UTC, is preserved in full. It reports: proper local-monoid sizes1/2 and the full local15; power profile12 fixed, two index-two collapses and one period-two unit; a14-state idempotent-generated sector with three minimum triples; all15 elements regular; only units6/8 left/right cancellative. These are completed finite-table facts, not new probability laws. Its statement that congruences/endomorphisms/shortlex were not touched by this five-task batch does not undo the merged #488–490 results.

### #484 head fixed once at705819e; N900 remains artifact-only here

The all-open list saw intermediate head `d23e3f94d0d5c03560d9ec77a8a379a1e47b18d7`. The final bounded #484 read is open/non-draft at **`705819e95d1146fdedb06e9c7628344f108b80af`**, with61 changed paths and still no discussion comments. This report does not chase later heads.

The new path is `notes/p334-general-two-port-birth-theorem.md`: d23e3f94 first adds190 lines under “Derive general two-terminal birth theorem and parallel reliability factorization”;705819e adds32 lines under “Derive exact winning-channel law from reliability polynomial derivatives”. These commit/path facts are captured. The note's scientific content is handed to the main reader for full review, not certified here from its commit subjects.

N900 status remains **artifact_only / runtime_unknown**: manifest and runner are published, but no N900 raw, result or receipt paths appear in the final captured PR changes. The body's older “not started” sentence is not a process observation. No process/server inspection or launch was performed.

### Active scientific discussion that navigation must consume

[Draft267 comment5474385016](https://github.com/LightChainr/Matching-One/pull/267#issuecomment-5474385016),06:01:15 UTC, has substantive completed science beyond the current Issue bodies:

- **P334 full conditional law is already solved**, at `6358ba49ef390c10a3f501b589ba7ba1d4e05b09`: two fixed N425 prefixes, all2^173 insertion subsets represented by two-terminal site connectivity, reported treewidth upper bounds4/6. Exact mean waits17.73237780/20.77877866; final A/B gap3.04640085; selected minimal quintics1141/7196. The same witnesses' site attribution is also complete at `1c06230b8f7e13be98f128361ad72b23c0c425ae`: B's leading-five-port share43.0173% falls to32.3775% in T>40, while internal-site share rises23.6534%→40.8146%. Thus “first longer-horizon/full reliability” must not remain a next task. The #492 dual-cycle carrier does not itself supply these laws; these are separately delivered results on the same dependency block.
- **N400 production, finite A/E secant and scale comparison are done**, as above. The same comment separately names the ordinary scalar/no-Jacobian comparison at `8ac86fb`/`162fe94e`: N1009.648/2 versus N400.07698/2, explicitly target-underpowered. That model is not the finite-density-coordinate model and is not a second vote on it.
- **P398 hidden-geometry memory is already analyzed**, at `39e06607ec3a353b1130acebf770da591acaf340` and `c9dc218f5522502cce8cca539b876ed5faa49b8a`. The comment reports R/T2 hidden-force support, memory times≈.272/.271, unfitted compression crossing.25412 against full.26566;6.812-fold normalized feedback is not bare coupling, since plus bare feedback is.5461 times minus while its source variance is≈12.47 times smaller. R/T2 are degenerate within each protected ray; the named cluster-contact hierarchy is not two independent fields or new Monte Carlo.

These are authored source reports preserved without rerunning them. They identify completed outputs that a current scientific map must consume; neither a source comment nor a title substitutes for the main agent's pinned-report review.

The limited active-Issue check found concrete stale current-body wording:

| Item | Current-window discussion | Correction needed in the scientific navigation |
|---|---|---|
| #154 / #275 | No comments after05:45:09 UTC; current bodies still accurately record c0880c2's completed but unresolved mixed row. | No newly completed distinct energy-field identification found in this increment; do not reassign the first local interaction row. |
| #370 | Issue, not PR; last updated05:09:38 UTC, no current-window comments. | Keep f5779b91's existing-production certificate separate from already-built exact/compiler support. Additional support work is not permission to analyze E_top. |
| #205 | No new comments;06:04 body says N400 manifest is not acquired data and proposes first scale transport. | Superseded by actual PR484 N400 data, finite secant and full-scale reports. |
| #334 / #429 | #429 adds comment5474342657 linking #492; #334 adds no comment. Bodies still propose longer-horizon reliability. | Preserve #492/#491 distinctions, but consume6358ba49/1c06230b before assigning new work: full conditional law and site attribution are already done. Population loading remains a different unmeasured object. |
| #398 | No new comments;06:04 body stops at fixed-readout/width8. | Later PR267 comment already supplies projected memory and contact-geometry decomposition; “first hidden-force interpretation” is stale. |
| #418 | No new comments;06:04 current body records corrected joint compatibility and unreliable r5/sharing. | No new scientific increment found; do not restore the old radius-flow claim or reassign normalized replay. |

This is an increment check of the active items returned since05:45:09 UTC, not a second all-Issue review. No issue was closed, locked, relabeled or otherwise changed.

### Existing open PR heads

The read-only open-PR listing returned23 entries (one Draft). Original titles and full heads are preserved below; except for #484's separately recorded final705819e overlay, this is the14:16 listing, not an assertion of future state. All target main except #277, whose base remains `results/p154-norm4-production-reveal-20260829`. No open-PR bodies outside this bounded task were reread just because they appear here.

| PR | Original title | Status | Head at open-list capture |
|---|---|---|---|
| [#492](https://github.com/LightChainr/Matching-One/pull/492) | P429: certify dual-cycle blockers on the two real N425 trigger graphs | open | `0e52dbaeed53dfffa94592e53e38129c179c5078` |
| [#491](https://github.com/LightChainr/Matching-One/pull/491) | Prove the rank-one cut-network theorem and explain the real trigger witnesses | open | `ab90201e88409310632812727e0138c56b455644` |
| [#485](https://github.com/LightChainr/Matching-One/pull/485) | N100 exact shape nulls and clock-orthogonal thermal deformation | open | `b454bb8ec04ad90a2db579efaf6285e59d6ba5aa` |
| [#484](https://github.com/LightChainr/Matching-One/pull/484) | N100 production: shape splitting and zero-area thermal redistribution | open | `d23e3f94d0d5c03560d9ec77a8a379a1e47b18d7` |
| [#451](https://github.com/LightChainr/Matching-One/pull/451) | Score real same-stream P439 crosswalk and unresolved M loading | open | `bfbceb24f4072e5fd2025a2cecb344014adbd9d8` |
| [#438](https://github.com/LightChainr/Matching-One/pull/438) | Exact: close the W5 relative-dual state and periodic gluing (#14) | open | `bd2561abd889b7360178616467422e31ad91c838` |
| [#435](https://github.com/LightChainr/Matching-One/pull/435) | Exact: full survival laws do not close the topology state (#429) | open | `ffd91ebd819f7893cbee84aeb3f40da14a700a7b` |
| [#417](https://github.com/LightChainr/Matching-One/pull/417) | Exact controls: generator-dependent rank, cooperative birth survival, and typed correction design | open | `5e047117de37724c8aa1cb56b3818a8018def13f` |
| [#416](https://github.com/LightChainr/Matching-One/pull/416) | Exact: P250 positive spatial spectrum and a ≥100-mode witness | open | `e105b00636a07cfe5c650f80ffe34f8cfc10ba81` |
| [#415](https://github.com/LightChainr/Matching-One/pull/415) | Certify cooperative homology continuation beyond one-step H2 | open | `d09f9252dd18de497ab88bd30e92101813432a81` |
| [#385](https://github.com/LightChainr/Matching-One/pull/385) | Exact identifiability controls: Jordan closure and a C4-protected five-state quotient | open | `b5c761a05019f33b5528ecb9c4aa6a18006e71c1` |
| [#277](https://github.com/LightChainr/Matching-One/pull/277) | Fit minimal Jordan plus one even-mode transfer after norm-4 reveal | open | `3e855ced4fd98d8979c0b712636b45c2fa54f969` |
| [#273](https://github.com/LightChainr/Matching-One/pull/273) | Reveal frozen norm-4 q=2 versus Jordan production score | open | `8b26a30a785bc142a9d17bfed99a8d0e98ddc4dc` |
| [#267](https://github.com/LightChainr/Matching-One/pull/267) | docs: recover scientific frontier and score production mechanisms | Draft | `7cdf1d575f077cebec5cdfb85a25695cc03d1728` |
| [#247](https://github.com/LightChainr/Matching-One/pull/247) | Run the norm-5 two-cutoff pivotal score | open | `9d35fa65de5bd4c81321751024dc8007dc329a87` |
| [#246](https://github.com/LightChainr/Matching-One/pull/246) | Freeze triangular energy/log-pair sufficient statistics | open | `7f9dcd88eb14cf89dce499902953ade1c57abcb7` |
| [#245](https://github.com/LightChainr/Matching-One/pull/245) | Add exact Boolean noise-semigroup oracles | open | `c7eddeb53879182287de4c555fad4312249f5ee9` |
| [#230](https://github.com/LightChainr/Matching-One/pull/230) | Freeze multi-u thermal-response templates from the #101 coordinate map | open | `9858eac23e3c8b568c8661684f22d5ea02704fc8` |
| [#229](https://github.com/LightChainr/Matching-One/pull/229) | Identify the matching polynomial as a vertex-subset defect generating function | open | `83a14f3d74295fc9218fdd5be4905980dd40742e` |
| [#228](https://github.com/LightChainr/Matching-One/pull/228) | Certify Galois groups of committed axis matching polynomials | open | `9ab422978e380d96516760014582ad7bf17a02b5` |
| [#197](https://github.com/LightChainr/Matching-One/pull/197) | Reanalyze local matching zeros near the physical root | open | `bf31f18429f2ef929fdc385e4efe0a2e4a726e42` |
| [#196](https://github.com/LightChainr/Matching-One/pull/196) | Lift matching identity to configuration Euler/Betti observables | open | `4ebfc91753e9569983b0a40b8d3eeefdeac0b172` |
| [#84](https://github.com/LightChainr/Matching-One/pull/84) | Research: exact axis L=5 matching polynomial frontier | open | `6a049b764658bc15956ae83f704c85e1d36ff660` |

**Final handoff:** this two-file capture is complete and released. The currently missing physical question is not “whether E_top has ever been tested”; it is how the still-unidentified original norm-4 response relates to the now measured topology/clock/local-source and finite-geometry structures. New finite support results can help that comparison without becoming a new preparation queue. No new experiment is launched by this report.

