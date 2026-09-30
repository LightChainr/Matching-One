# Exact averaging resolves a large-size completion-creation history contrast

2026-09-30. Same 140,000 square-L512 prefixes as the
[safe-insertion block](safe-insertion-independent-block-20260929.md), not a
new independent sample. The target, entry cutoff, exact (D,c) cells and
overlap weights are unchanged. The estimator was improved after the first
readout; this timing and the original result are retained.

**Outcome:** early births create fewer new completion sites than late births
at the same current direction and completion count. Exact conditional
averaging gives **−0.00344939 +/- 0.00077225** (aligned 14-batch SE), versus
the original noisy **+0.00668183 +/- 0.01208706**. The uncertainty shrinks
by a factor **15.65**, revealing a finite-size history signal about 4.47 SE
from zero. This is substantive large-size evidence against recursive
closure of the enlarged (rank,D,c) observer, not a proof of its continuum
behaviour or a second independent confirmation.

## 1. What changed computationally

For each current occupied set A, the new
[projected-component criterion](projected-component-pair-criterion-20260929.md)
enumerates the exact unordered safe synergy pairs. Shared occupied
components impose integer transverse-potential constraints; disagreements
between them, or with a direct physical edge, detect exactly the pairs
that complete together. The fixed-degree expected cost is O(N+e), where e
is the number of output pairs, not an all-vacancy-pair scan.

With m=N-b, c=nu_b and e pairs, the mean increment of a uniform safe
insertion is exactly 2e/(m-c). Replacing four random probes with this
conditional expectation removes probe randomness. The original prefix
sample, all cell weights and finite-prefix uncertainty remain. This is
conditional averaging of the same target, not oversampling high-contact
sites or changing the preparation source.

All 140,000 tuples `(J1,rank_b,dx_b,dy_b,nu_b)` matched the original data
row by row. Integer c/e values are saved. Reported cohort estimates and
covariances are floating-point statistical calculations, not exact
rational population probabilities.

## 2. Actual complete-block result

Early means J1<=154646; late means 154646<J1<=155385; rank_b=1,
L=512, N=262144. Exact-cell weights are n_E*n_L/(n_E+n_L).

| Early minus late, unless otherwise marked | Estimate | Batch SE |
|---|---:|---:|
| Exact conditional mean safe increment | −0.00344939 | 0.00077225 |
| Original four-probe increment | +0.00668183 | 0.01208706 |
| Paired exact-minus-probe estimator difference | −0.01013121 | 0.01239566 |
| Exact mean synergy-pair count contrast | −183.96398 | 41.20529 |
| Derived two-next-insertion survival contrast | +3.22818e−8 | 7.23066e−9 |
| Exact weighted early increment mean | 0.09557456 | 0.00045022 |
| Exact weighted late increment mean | 0.09902395 | 0.00052572 |
| Unchanged secondary: same-D current c contrast | −0.03986512 | 0.56987795 |

The mean increment difference is about 3.48% of the weighted late mean;
this is a relative response, not an explained or mediated fraction.
The original and refined point estimates differ by only about 0.82 SE
of their paired difference: the old positive point was poorly resolved,
not an established opposite physical effect.

The original primary and secondary points and SEs are reproduced by the
aligned data, without rerunning their original scoring pipeline. This is
an identity check, not independent scientific replication. Every deletion
removes the same entire batch from both measurements and recomputes support
and weights. The result retains the full 8x8 covariance and 14 deletion
vectors, including the exact dependence among the derived quantities.

There are 55,028 rank-one prefixes. The same 1,190 common (D,c) cells retain
54,463 of them, 98.973%; no binning or replacement definition was introduced.
The empirical SE ratio is 0.063891 and variance ratio 0.004082. Those numbers
describe this estimator comparison; they are not an increase in independent
sample count or a claim that every future jackknife must improve equally.

## 3. Scientific interpretation and next distinction

At fixed c the immediate exit hazard is identical by definition. The
negative edge-count contrast gives a positive two-insertion survival
contrast by the exact transfer −2*Delta(mean e)/[m(m-1)]. Thus this signal
has a concrete future-rank consumer; it is not just a geometric feature
correlated with age. A nonzero population contrast would violate a
necessary condition for weak Markov closure of (rank,D,c), because the
birth cohort is part of that observer's past. This finite-sample result
supports such a violation at L512 under the stated target.

It does **not** revive the older direction-only claim: its independent
secondary remains −0.040 +/- .570 and did not reproduce +1.629 +/- .391.
Adding c changes the observer. The enlarged observer can expose memory
even when a coarser observed process is Markov; these two statements are
not contradictory. Nor does one nonzero mean identify the full successor
law, the minimal state, original U, or a continuum operator.

Next attention goes to two concrete questions, not another feature scan:

- Fix this same negative-sign completion-creation prediction and score it
  once on an independent prefix block with the exact evaluator. Prefer a
  correctly replayable independent archive before collecting new samples.
  If it is another previously inspected archive, say so; do not call it
  prospective new acquisition. Do not reuse the present block as the second
  confirmation or pool it to manufacture a stronger number.
- Derive a continuous-kernel prediction for creation versus survivor
  selection. The [new label-time identity](completion-pair-label-curvature-20260930.md)
  gives hazard'=[mean c+2 mean e−Var(c)]/(1-t)^2. Fixed-count data cannot
  be substituted by t=b/N, and a two-insertion effect cannot be extended
  to a 735-insertion window by assuming constant curvature. Scaling and
  survival-selection cancellation are the actual next mechanism questions.

The finite L4 counterexamples remain exact; this is separate larger-size
evidence, not an asymptotic extension of their proof. These are priorities,
not locks on parallel research.

## 4. Execution and reproducibility

Eight completed batches were retained locally. After the owner's request
to prioritize Huawei for medium/large work, local incomplete workers were
stopped and the remaining six batches ran on TV2N0X, six workers, in
420.381 seconds including compilation and a 16-old-prefix cost check.
Actual available budget was 14.5 CPU / 25 GiB; GCC10.3.1, Python3.9.9,
ARM64. Other machines were not operated. The prior user pause and abandoned
partial runs are not included in successful-batch runtimes and do not add
samples. Their completed checkpoints were preserved.

The engine/geometry code was unchanged between local and cloud evaluation.
Replay source: `8b5eee0d47b6c220c2ba46a5bae08d2d19b392a4`; cloud handoff:
`0684c3b4eee0a647e643540ecfcdcd2918a70977`. Four small physical pair-set
controls and a 128-old-prefix identity/cost check ran once before full
replay. The scorer's small two-cell arithmetic control and one complete
production scoring pass were performed; no full suite or old census rerun.

[Run and raw rows](../analysis/exact-completion-pairs-20260929/data/run.json) ·
[Results and covariance](../analysis/exact-completion-pairs-20260929/results/result.json) ·
[Compact table](../analysis/exact-completion-pairs-20260929/results/RESULT.md) ·
[Method and reproduction](../analysis/exact-completion-pairs-20260929/README.md).
The old four-probe report is preserved unchanged.
