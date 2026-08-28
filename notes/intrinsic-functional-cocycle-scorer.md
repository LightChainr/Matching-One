# Intrinsic functional-cocycle scorer

This implements the pre-target full-curve discriminator proposed in Issues #101, #119, #125 and #138. It is deliberately narrower than a new scaling fit: the unknown leading function and the unknown first correction function are eliminated algebraically at each already-frozen intrinsic level `u={0,.025,.05}`.

## Frozen computation and canonical entrypoint

The 2026-08-28 kernel and scoring contract are retained unchanged:

```text
scripts/score_intrinsic_functional_cocycle.py
predictions/intrinsic_functional_cocycle_scorer_20260828.yaml
```

After the type-safe channel layer in PR #170, the canonical target entrypoint is

```text
scripts/score_intrinsic_functional_cocycle_typed.py
```

with the additional pre-target semantic contract

```text
predictions/intrinsic_functional_cocycle_semantic_gate_20260829.yaml.
```

This adds no model and does not change q=2/Jordan chronology. It only makes the primitive full-curve observable semantics executable and fail-closed.

## Frozen functional score

For each size, the kernel solves `Mbar(p)=+/-u` inside every delete-one-batch replicate, reconstructs

```text
T_N(u) = D_even(u,N) + S_odd(u,N),
Z_N(u) = N^(13/8) T_N(u),
```

and scores

```text
R_c = Z_(5N) - c Z_(2N) + (c-1) Z_N
```

for `c=8/5` (ordinary relative-q=2 correction) before `c=log(5)/log(2)` (rank-2 Jordan cocycle).

## Primitive observable semantics

All source and target curves in this score come from the threshold-rank **cross** channel. The semantic gate freezes the primitive descriptors:

```text
Mbar center:  cross / odd  / raw value
P4[S]:        cross / even / angular-normalized contrast
P4[D]:        cross / odd  / angular-normalized contrast
P4[S']:       cross / even / angular-normalized contrast
P4[D']:       cross / odd  / angular-normalized contrast
```

Every P4 is divided by its own signed `DeltaCos4` at the same N before any cross-size comparison. Therefore there is no `either -> cross` conversion and no raw angular sign transported between sizes. Each primitive descriptor maps identically to itself with scale `+1`, offset `0`.

The typed wrapper validates all identity maps before calling the frozen kernel and appends the descriptor set plus applied maps to the output JSON.

## Covariance boundary

Existing N65/85/130/170 curves share counters and form one synchronized jackknife group. The N325 and N425 production runs use disjoint counters, so each is an independent group. The scorer sums those independent-group covariance contributions rather than manufacturing covariance from coincident numeric batch labels.

Every intrinsic center, `p_+/-`, `Mbar'`, P4 projector, center diagnostic and functional residual is recomputed inside the appropriate delete-one replicate.

## Production invocation

```bash
python3 scripts/score_intrinsic_functional_cocycle_typed.py \
  --histograms n65.hist.csv n85.hist.csv n130.hist.csv n170.hist.csv \
               n325.hist.csv n425.hist.csv \
  --covariance-groups 65,85,130,170 325 425 \
  --json functional-cocycle-score.json
```

The same run also reports the intrinsic-center diagnostics `J=P4[S']/Mbar'`, `N^(13/8)J`, and `Xi=J/P4[D]` with joint covariance. Those diagnostics test whether the bare-p thermal metric explains the center S-prime drift; they do not alter the frozen functional-score order.

## Claim boundary

The semantic gate does not change:

- the frozen u grid;
- either Gaussian lineage;
- the q=2 then Jordan model order;
- `8/5` or `log(5)/log(2)`;
- covariance grouping;
- the delete-one recomputation rule;
- the pseudoinverse cutoff;
- target data or fitted parameter count.

It only makes the full-curve source/target channel and normalization semantics explicit before the N325/N425 target score.
