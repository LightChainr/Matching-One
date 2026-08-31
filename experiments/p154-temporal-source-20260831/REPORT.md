# Earlier cluster geometry predicts later rank selection

**Completed on the original 2.4M permutations, with no new counter.** A source centered within its **earlier** rank has zero direct effect on that earlier rank distribution, yet it strongly changes subsequent first and second activation probabilities on all six norm4 tori. Thus early rank alone is not a sufficient predictive state for this specified cluster-count source and lag. The original H4 readout U still has no resolved source response.

The new information is temporal: `s_L` is observed at `L=max(0,K-ceil(sqrt(N)))`, while q/E and entry/exit are observed at later K. The earlier equal-time source, line, soft-angle and occupancy/rank projections did not contain these products. This calculation freezes PR267 `7da1eeb0e51cf430987dbf204d23713c2ab5a46c` and original engines `bfab0330f5f56ca4d746b45d737f1607e3d229a0`. It does not repeat their primary analyses.

## A finite mechanism prediction that fails

At fixed L and geometry, take `s_L=C_B+C_W` and

```text
z_K=s_L-E[s_L | L,r_L,g].
```

For any function of the early rank, `Cov(f(r_L),z_K|K,g)=0`. A proposed rank-only predictive state would additionally imply

```text
E[f(r_K)|early configuration] = E[f(r_K)|L,r_L,g]
    => Cov(f(r_K),z_K|K,g)=0.
```

The measured later response violates that implication. Conditional cluster geometry retained by s_L carries information about future activation beyond early rank. Under the declared positive path tilt `exp(t*z_K)`, more fragmented early configurations reduce later cumulative activation probabilities at the root window. This is a finite ensemble response, not a fitted transition-rate law.

The immediate invariance is a **direct source derivative at fixed p/K**. Following the moving matching root also changes the occupancy mixture; no claim is made that the early-rank population remains invariant under that additional thermal transport.

## Later activation has two distinguishable contributions

All quantities below differentiate the bulk source strength t. Errors are one paired-jackknife standard error. At N260 (lag17), the two geometry responses are:

| Centered early source stratum | First entry, first geometry | First entry, second geometry | Second exit, first geometry | Second exit, second geometry |
|---|---:|---:|---:|---:|
| Early rank0 | −.228771±.000879 | −.229618±.000917 | −.165941±.000670 | −.167796±.000663 |
| Early rank1 | exactly0 | exactly0 | −.047207±.000350 | −.047466±.000397 |
| Early rank2 | exactly0 | exactly0 | exactly0 | exactly0 |
| Sum | −.228771±.000879 | −.229618±.000917 | −.213148±.000774 | −.215262±.000831 |

At N340 (lag19), total first-entry responses are `−.248739±.001007` and `−.250348±.000958`; second-exit responses are `−.233835±.000891` and `−.236005±.000855`. All twelve geometry-specific centered entry responses and all twelve exit responses are negative and resolved by at least70 SE across the six sizes.

The early-rank1 contribution can only affect exit: first entry has already occurred. Its exit response is negative at every size and geometry, resolved by at least29 SE. Early rank2 is absorbing. These structural zeros are exact consequences of monotonic rank and centering, not independent null experiments.

Consequently, following the root, the early-rank0 part decreases rank1 population while the early-rank1 part increases it. At N260 their contributions are `−.06240388±.00047696` and `+.04732701±.00027333`, summing to `−.01507687±.00059341`. At N340 they are `−.07145602±.00064900` and `+.05683695±.00028717`, summing to `−.01461907±.00077046`. These are correlated components; their errors must not be combined as independent errors. Their partial cancellation is measured rather than assumed.

Entry/exit here mean cumulative events `I[K>=K_minus]` and `I[K>=K_plus]`. They are not event-density scores, hazards or estimates of a universal stochastic ordering at every K.

## Original U and its complete moving-root derivative

The original U uses the same exponent13/8, exact geometry contrast, pooled matching root and common slope as the earlier norm4 analysis. Its derivative includes all four terms: direct response, root motion, source-dependent slope and slope transport. The centered early-source results are:

| N | Fixed lag | Root derivative | Root-comoving rank1 derivative | Original U derivative v |
|---:|---:|---:|---:|---:|
| 65 | 9 | +.0250181±.0001925 | −.0170792±.0009968 | −.13805±.85203 |
| 85 | 10 | +.0271196±.0001661 | −.0163782±.0011606 | −.04607±1.23867 |
| 130 | 12 | +.0292956±.0002283 | −.0155022±.0015128 | −1.94282±4.33670 |
| 170 | 14 | +.0308001±.0002299 | −.0151191±.0016087 | −3.04013±5.98652 |
| 260 | 17 | +.0315499±.0000776 | −.0150769±.0005934 | +.84255±4.88218 |
| 340 | 19 | +.0312034±.0000675 | −.0146191±.0007705 | +12.24934±9.92174 |

Root motion and root-comoving population response are resolved in all six cases. Every centered v is within1.24 SE of zero. The q2 and Jordan contrasts of v remain unresolved on both lineages: q2 is `7.37550±16.30145` / `33.57299±26.83415`; Jordan is `4.59014±10.00300` / `18.28353±15.64654`. These are source-extension diagnostics, not new tests of the already rejected high-statistics unperturbed q2 model.

The uncentered early source is also fully reported in [latest.json](results/latest.json). It gives endpoint v `−1.07994±5.05984` and `14.83646±11.02107`, likewise unresolved. No weak total is converted into a percentage allocation or used to select a new continuum field.

## What was acquired and checked

Exactly100k old permutations at N65/85/130/170 and1M each at N260/340 were reobserved once, with the two geometries sharing each permutation. The endpoint batches retain their old1000 plus additional9000 union. Every original batch/K q,E,s,qs,Es sum was reconstructed **exactly as integers** before scoring. This verifies the counter mapping, geometry conventions, source units and batch alignment alongside the new observations.

The final670-dimensional joint covariance retains the original three source groups and100 aligned deletions. Every deletion re-estimates the early conditional means and reuses the correct saved pooled root. Central and omitted old equal-time controls agree with the previous results to `2.12e−12` for U and `6.32e−11` for v. Centered stratum addback error is at most `2.81e−12`; early-rank1 first-entry source sum is exactly0. Other theoretical zeros have raw-sum rounding residual at most `3.73e−9` before division by100k/1M.

A separate [numerical derivative check](results/moving_root_verification.json) solves roots of the auxiliary first-order `(1+t*z)` measure at symmetric t steps for both sources and all six N. Its derivative differs from the four-term formula by at most `1.83e−7`. Those48 small root solves validate the derivative; they are neither new simulations nor finite-strength experimental results.

The assigned ZyTrST machine used16 threads, GCC10.3.1, Python3.9.9, NumPy1.26.4 and SciPy1.13.1. Acquisition including compilation, exact checks and compression took41.644s; the acquisition-plus-analysis driver finished in50.200s. The source score itself took7.249s. [run.json](results/run.json) records each command, duration and hash. Results were retrieved and their hashes verified before requesting VM shutdown. The final lifecycle record is in [RUN_RECORD.md](RUN_RECORD.md).

## Scope and next scientific use

This is a **new lagged path-source perturbation**, with its earlier observation time tied to K. It is not a decomposition of the old equal-time source derivative. The centering coefficients are empirical, geometry-specific and re-estimated; they provide a conditional predictive mechanism test, not an independently frozen physical intervention. In p differentiation, the lag, K-to-L mapping and conditional means are held fixed while only Binomial weights vary. No source-time derivative or new exponent is introduced.

The result rules out the stated early-rank-only predictive sufficiency for this cluster source and selected lag, and identifies a resolved early-rank1-to-exit channel. It does not identify the original norm4 field, a full transition kernel, an asymptotic temporal exponent or a finite-strength source map. The same archived permutations underpin all outputs.

The next mechanism should explain or predict the measured **entry suppression versus exit suppression**, especially the opposing early-rank0/rank1 population contributions, and then predict the original geometry contrast. A wider lag sweep or another first signal-existence run is not required to establish the channel demonstrated here.
