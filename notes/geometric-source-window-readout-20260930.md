# L512: an invisible immediate source has a resolved735-step rank response

2026-09-30. The [fixed finite-window contract](geometric-source-window-identity-20260930.md)
has now been evaluated, not left as another tool or production proposal.

**Result.** A single clock-preserving geometric kick has a clearly resolved
response after the original735-insertion window. It does not disappear
between local completion creation and this endpoint. The source preserves
the current exit probability, so this is an actual delayed geometric
transmission, not a change of the occupation clock.

| Conditional on original rank_b=1 | Estimate | Aligned14-batch SE |
|---|---:|---:|
| d/dtheta survival at b+735, theta=0 | **-0.01104203191** | **0.00015766645** |
| Baseline survival at b+735 | 0.4447553973 | 0.0021409587 |
| Exact two-step source susceptibility | -0.000059071990 | 0.000000469939 |

The735-step derivative is about70 batch SE from zero. These are derivatives
per unit theta in exp(theta*degree), not finite-theta probability changes.
All three coordinates are correlated readouts of the same original block;
the full3x3 covariance and14 deletion vectors are retained. The two-step
quantity is a geometric identity, not a second independent experiment.

## What was actually evaluated

- Square NN occupied-site torus, L512, N262144; unchanged b155385 and
  endpoint156120. Exactly one735-step window, no lag/angle/size/feature scan.
- The original safe-insertion block's140000 complete Fisher-Yates
  permutations were replayed. There are55028 rank-one prefixes. All140000
  first-five-field tuples match the earlier raw prefixes.
- At each rank-one prefix A, enumerate its exact synergy graph once,
  then follow the **original uniform continuation** through735 arrivals.
  The estimator is the exact conditional score identity

      1{rank_{b+735}=1} [sum_{v in next735} degree_A(v)/735 - 2e_A/s_A].

  It averages over which of that unordered arrival set came first. It
  does not recompute initial degrees later, simulate a finite-theta model,
  or extrapolate the two-step curvature.
- Conditional risk denominators are recomputed on each whole-batch
  deletion. No early/late subgroup was added to this new contract, and no
  precision or sample extension was triggered by the result.

This is a new mechanism response on **previously analysed prefixes**, not
independent evidence for their passive birth-history contrast. The separate
[independent archive reproduction](completion-creation-independent-replication-20260930.md)
is already recorded and is not pooled with this calculation.

## What now follows, and what still does not

1. **The finite-window channel exists at this large size.** The selected
   micro-source is invisible to immediate rank and to the occupation clock,
   but its effect reaches future rank after735 insertions. A pure clock
   explanation of this intervention is insufficient.
2. **This is a second-birth actuator.** Since the kick acts only after
   rank one is reached, the first-birth distribution is unchanged. Thus at
   the fixed count endpoint, unconditional deltaM=deltaE_top=deltaP2, with
   conditional survival response multiplied by minus the original rank-one
   risk probability. It is not an observed change of the first-birth source.
3. **The result does not yet attribute the natural memory signal.** A strong
   response to an intervention does not determine how much of the passive
   early/late association arises at birth, during safe transport, or through
   survivor selection. That separation is the next mechanism calculation,
   using the already derived safe-path/weight identity, not more precision
   on this response.
4. **No label-time or original-U substitution.** A count endpoint is not a
   fixed p. The [one-sided source map](one-sided-birth-source-normal-response-20260930.md)
   proves what survives moving-root projection and supplies an exact L4
   symmetric positive control; a numerical L512 label-time map still needs
   its proper count mixture. No continuum exponent, universal amplitude,
   physical H4 attribution or original norm-4/U identification is claimed.

In particular, this closes “does the supplied source transmit at all beyond
two steps?” for the fixed L512735-step experiment. Do not automatically
repeat it at higher precision or start a size ladder. The useful next
question is source/selection attribution or the appropriate endpoint
mixture for a specified consumer.

## Execution and artifacts

TgFr7R (`83750ac4eae34bbbbc64b894b854dd8c`) executed14workers in519.145seconds,
with ARM64/Python3.9.9/GCC10.3.1, live14.5CPU and25GiB limits. A32-old-prefix
cost/correspondence check carried zero additional scientific weight.
The compiler was installed in the task's own installroot because the
restarted base image lacked g++; other workspace directories were untouched.

One support worker reviewed the exact score/source correspondence and one
specified L4,h3 rational control; the main investigator did not repeat
that verification. After transfer, one standard-library score in the local
research Python3.11 environment computed the result and all-row alignment.
No full suite, old census, fresh stochastic block or numerical theta grid.

[Result and covariance](../analysis/geometric-source-window-20260930/results/result.json),
[brief report](../analysis/geometric-source-window-20260930/results/RESULT.md),
[execution](../analysis/geometric-source-window-20260930/execution.json),
[source review](../analysis/geometric-source-window-20260930/verification.md).
Raw compressed rows, batch commands, exit statuses and source hashes are
retained in the same directory. Final resource release is recorded in the
execution file, not inferred from an ended SSH connection.
