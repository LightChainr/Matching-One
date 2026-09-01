# A spatial source one step earlier changes the next topological activation

## Technical summary

The new one-activation-lag experiment resolves a forward spatial mechanism
in the original norm4 square-site geometries. A cluster-count source centered
within the **previous occupancy and rank** changes the next activation even
though its same-time global response is zero. Root-comoving rank1 probability
decreases at all six N, by .00882 down to .00297 per unit bulk source strength,
resolved at15.5–37.8 batch standard errors.

Both activation channels contribute positively to the pooled matching-root
shift. The first-activation contribution is larger, and the opposing entry
and exit contributions leave a smaller negative rank1-population response.
This is a source-sensitive transition mechanism beyond the declared early
rank-only response model, not a second continuum-field identification.

The original-U H4 source response remains unresolved: at N260/N340 it is
**+.37897±.72219 / +1.27653±1.18227**. The mechanism reaches future rank
population; a directional mechanism for original U still needs a prediction.

One115.28-second local replay added joint marks to **2.4M old permutations**,
with zero new random samples. A1.364-second aggregation retained the original
three dependency groups and a96-dimensional joint covariance. Existing
matching curves were reconstructed exactly batch by batch; no root finder
or research test suite was rerun.

## The protocol changes the source time, not the baseline curve

For a fixed final occupancy K≥1, let j=K−1. In geometry g define the bulk
cluster count `s_j=CB(A_j)+CW(A_j)` and

```text
epsilon_j = s_j − E[s_j | j, rank_j, g].
```

Keep the early j/rank probabilities unchanged, exponentially tilt the
configurations within each layer, then add one uniformly chosen vacant site.
Read q and E at K, with `q=rank−1`, `E=q²`. Finally mix K using the original
`Bin(N,p)` weights; define the K=0 tangent to be zero. At zero source this
returns exactly the original uniform K-configuration law and matching curve.

This is a **one-step-before-readout injection** for each K. It is not one fixed
early injection shared across every later time, and not the original equilibrium
bulk fugacity acting on the final configuration. The full protocol, equations
and source units are in the
[pre-computation design](../../notes/norm4-lagged-source-mechanism-20260831.md).

Early conditional centering forces `E[epsilon_j q_j]=E[epsilon_j E_j]=0`.
It does not force the corresponding future moments to vanish. The new data
contain s at the **same permutation's K−1 prefix** paired with its activation
event; no product of previously saved marginal averages is substituted.

## Entry and exit oppose in rank population but reinforce the root shift

The six root-comoving rank1 responses, ordered N65,85,130,170,260,340, are

```text
−.0088170 ± .0002924     −.0068719 ± .0002933
−.0054583 ± .0003146     −.0046278 ± .0002979
−.0035436 ± .0000937     −.0029708 ± .0000919
```

Their standard-error ratios are −30.15,−23.43,−17.35,−15.53,−37.80,−32.33.
The four smaller N share a random-stream group; these are not six independent
replications or six new pieces of evidence to multiply.

At N260, the pooled root shift is `+.002773042±.000006763`.
Its first-activation contribution is `+.001512320±.000004709`, and its
second-completion contribution is `+.001260722±.000004760`.
At N340, the corresponding values are
`+.002321931±.000005780`, `+.001256648±.000003885`, and
`+.001065284±.000004366`. These are additive contributions to the **same
matching-root displacement**, not separate median birth-time shifts or a
claim that every part of each birth distribution moves later.

The population response exposes competing effects. At N260, the exclusive
`0→1` contribution is `−.02098711±.00006537`, whereas `1→2` contributes
`+.01744361±.00006693`; their correlated sum, including the tiny direct-event
root-transport term, is `−.00354361±.00009374`. N340 similarly combines
`−.01930883±.00006044` and `+.01633798±.00006859`.
The net response is smaller than either resolved contribution. The shared
covariance is essential to interpreting the cancellation.

## The three event kernels are the mechanism-level output

For final K, let T01, T02 and T12 be the empirical means of epsilon at K−1
times the indicators of `0→1`, `0→2`, and `1→2`. All three use the full
retained permutation denominator. They give

```text
Hq(K) = T01(K) + 2 T02(K) + T12(K)
HE(K) = −T01(K) + T12(K)
HF1(K) = T01(K) + T02(K)
HF2(K) = T12(K) + T02(K).
```

Binomially mix these kernels to obtain the new source responses Jq, JE and
their analytic p derivatives. Root transport and slope normalization are
then recomputed for this source at the saved unperturbed roots. The old bulk
source's root derivative is not copied into this experiment.

Direct `0→2` events are present:5853/5854 of the two N260 million-permutation
archives, and4787/4641 at N340. They cancel from HE but contribute twice to
Hq, so dropping them would change the root shift and the U normalization.
Their rank1 effect here is only via root transport and is small; that is
not a reason to erase their q contribution.

There is a direct hazard interpretation. At fixed early j/rank, let h_ab(A)
be the fraction of currently vacant sites that would cause transition a→b.
Then T_ab is the early-rank probability times the conditional covariance
of s with h_ab. The measured response therefore identifies coupling between
within-rank spatial cluster structure and the next-step topological hazard.
A model whose specified source response is determined by early j/rank alone
would make these centered couplings zero. This statement does not require
the entire unperturbed rank process to be non-Markov in every sense.

## Directional U remains a separate target

Using the original exact DeltaCos4 and `N^(13/8)` normalization, the six
lagged-source U derivatives are

```text
N65    +.08902 ± .23984       N85    +.67558 ± .30909
N130   −.26292 ± .79113       N170  +1.58419 ± .95434
N260   +.37897 ± .72219       N340  +1.27653 ± 1.18227.
```

The N85 value is a2.19-SE hint in a correlated multi-readout exploration,
not a resolved common directional mechanism. No exponent was refitted, no
lag was scanned, and the original q2/Jordan extension hypotheses were not
automatically applied to this different source protocol. Strong population
response cannot replace the missing directional identification.

## Data, uncertainty and reproducibility

The source definitions and code were committed before marking and scoring at
`6ac1d1379414b5788d94f270d174e8611a54e959`.
The [contract](../../analysis/norm4_lagged_source_contract.json),
[replay driver](../../scripts/replay_norm4_lagged_source.py) and
[analysis](../../scripts/analyze_norm4_lagged_source.py) retain:

- Four cyclic N use the original shared100k counters; each endpoint uses
  its original1M counters and its own seed. Both physical orientations stay
  paired; N260/N340 retain the actual HNF quotient from the frozen engine.
- Each endpoint batch is its original1000 plus later9000 permutations,
  not a new consecutive10000 block. The raw replay is stored as six gzip
  CSV files totalling about1.23MB, with hashes in [run.json](run.json).
- Each aligned delete-one removes the same event and early-source rows,
  reestimates early-rank conditional means, and uses the matching retained
  sample's saved root. The JSON retains all96 readouts, full covariance and
  each group's100 delete-one vectors and covariance factor.
- New event counts reconstruct original q/E sums **exactly for every batch,
  geometry and K**. Additive response reconstruction is within2.08e−13.
  These checks concern this delivered calculation, not a repeated test suite.

The centering uses plug-in conditional means and propagates their uncertainty
through the same jackknife. Reported SEs are finite-sample exploratory estimates,
not exact confidence certificates. The zero same-time response is imposed by
centering; the future-event covariance is the new measurement. No independent
sample population, unique microscopic field or continuum limit is claimed.

Reproduce in a separate worktree at the execution commit, before the saved
outputs exist, with the pinned backend objects available:

```bash
python scripts/replay_norm4_lagged_source.py --workers 4
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python scripts/analyze_norm4_lagged_source.py
```

## Next prediction: how winding contacts organize the hazard

The first lagged-source measurement is complete. A productive next model
should predict the **signed T01/T12 amplitudes and their direction contrast**
from which winding components can meet at a vacant site. The existing kernels
provide the target; another same-time spatial detection or an unrestricted lag
scan would not specify that mechanism.

One concrete hypothesis is that the source preferentially selects fragmented
within-rank configurations whose winding-contact opportunities suppress entry
more than completion. The measured root and population signs motivate it;
the relevant contact/winding counts have not been measured here. It remains
a hypothesis until a declared contact model predicts those kernels or their
residual directional response. P334 contact loading and P398 incidence counts
are useful conceptual guides with distinct sources and state spaces, not
substitute evidence for this norm4 prediction.
