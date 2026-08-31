# P154: early cluster source and later rank selection

This package adds one lagged source/rank joint observation to the **same 2.4M old permutations** used by the completed six-size norm4 source analysis. It adds no random counter and no new geometry. The sole lag is fixed before acquisition: `ell=ceil(sqrt(N))` (9,10,12,14,17,19 in increasing N order).

The baseline code and results are frozen at PR267 `7da1eeb0e51cf430987dbf204d23713c2ab5a46c`; the exact production engines are at `bfab0330f5f56ca4d746b45d737f1607e3d229a0`. [SOURCES.json](SOURCES.json) identifies every imported object. The previously completed occupancy/rank projection at `8799dfe1` concerns an equal-time source. This experiment asks a different, lagged question with newly observed joint marks.

## Scientific object

For a permutation `omega`, let `r_j` be the black ambient rank after j occupied sites, and let `s_j=C_B(j)+C_W(N-j)` be its bulk black primal plus complementary white matching component count. At late observation occupancy K, set

```text
L(K)=max(0,K-ell),
z_K=s_L-E[s_L|L,r_L,g],
q_K=-1+I[K>=K_minus]+I[K>=K_plus],   E_K=q_K².
```

The centering is within early rank and geometry at fixed source occupancy, not within late rank. It ensures zero response of every function of `(L,r_L,g)` and every K-only observable, but imposes no zero-response identity on later q/E. The already completed equal-time projection cannot recover `E[q_K s_L]` or `E[E_K s_L]`; the old profiles contained only `q_K s_K` and `E_K s_K`.

The specified positive path measure is, in each geometry,

```text
P_(p,t)(K,omega) proportional to Bin(N,p)[K] * P_uniform(omega) * exp(t*z_K).
```

The uncentered early source `s_L` is retained as a control. All derivatives use bulk t; there is no extra factor N. When p changes, the lag, mapping `L(K)` and fixed-K conditional means remain fixed. Only the Binomial weights change. Thus `d_p` already includes the changing distribution of source and observation occupancies; differentiating the step function L or the fitted conditional means would define a different source.

Let `Jq=Cov(q_K,z_K)` and `JE=Cov(E_K,z_K)` under the above baseline. Since centered source has zero mean at every K,

```text
Jq = sum_K w_K sum_r [ E(q_K*s_L*I[r_L=r]|K)
                    - E(q_K*I[r_L=r]|K)*E(s_L|L,r,g) ],
```

with the identical formula for JE. Empty empirical strata contribute zero without imputation. `d_p Jq` and `d_p JE` differentiate weights only.

The original readout and its moving-root derivative are retained:

```text
D=mean_g q_p, B=P4[E_p], H=P4[E_pp], T=mean_g q_pp,
U=N^(13/8) B/(2D),                  p0dot=-mean_g Jq/D,
v=N^(13/8)/2 * [ (P4[(JE)_p]+p0dot*H)/D
                 - B*(mean_g (Jq)_p+p0dot*T)/D² ].
```

`P4` uses the unchanged first-minus-second geometry contrast and exact Delta cos4: `1152/845` or `2304/1445`. No exponent is refitted. Each paired jackknife uses its previously saved pooled root and recomputes the early conditional means, U jets, root motion and slope response.

## Named entry/exit channels

At fixed p, first entry and second exit respond as

```text
J_entry=(Jq-JE)/2,   J_exit=(Jq+JE)/2,
J_rank1=J_entry-J_exit=-JE.
```

These are cumulative first/second activation events at the later K, not event-density derivatives or a hazard model. The centered source is split into early-rank strata before scoring. Exact monotonic-rank consequences provide nontrivial semantic checks:

- early rank2 is absorbing, so its later q/E response is zero;
- in early rank1, first entry has already occurred, so only the later exit can respond;
- early rank0 can influence either later first entry or second exit;
- the three centered strata add to the centered source's v, root motion and rank1 population response, with their full shared covariance.

## Inputs and replay identity

`inputs/old_profiles.npz` losslessly preserves the integer sums `(q,E,s,q*s,E*s)` for every original batch, geometry and K; it is derived from the sources listed in SOURCES.json. `inputs/anchors.json` extracts each central/delete-one root and original same-time source readout from the frozen endpoint analysis. Vendor files are verbatim production engines, retaining original RNG, quotient labels and homology convention.

For N65/85/130/170, seed2026104501 uses original counters `[5100000000,5100100000)`. For N260/340, seeds2026105401/2026105402 use `[8200000000,8201000000)`. Every endpoint analysis batch b remains its old1000 counters at offset1000b **union** its additional9000 at offset100000+9000b. Both geometries use the same permutation. The four cyclic sizes preserve their common dependence group; each endpoint remains its own group.

`replay.cpp` also reobserves the already known equal-time five moments in the same two sweeps. `run.py` requires **exact integer agreement at every batch/K**, before accepting the new marks. This is a validation of the same acquisition, not a second replay. It writes compressed raw moments and execution hashes, then `analyze.py` propagates the original3 groups ×100 paired omissions.

## Reproduce

On Linux with C++17/OpenMP, Python, NumPy and SciPy:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=16 python3 run.py
```

Run in a fresh copy without an existing `results/run.json`; the driver refuses to repeat acquisition. To repeat only the analysis from existing raw marks, use a separate copy without `results/latest.json` and run `OPENBLAS_NUM_THREADS=1 python3 analyze.py`. Production results are never silently overwritten. The package needs no Git checkout or network access to run.

## Interpretation boundary

This is a new source defined on a path, with source occupancy tied to the observation occupancy. It is **not an allocation of the old equal-time physical source derivative**. Empirical conditional centering is re-estimated and descriptive; the calculation does not claim an independently fixed microscopic intervention or causal identification. Lagged entry/exit responses are finite predictive geometry channels under a positive tilted path measure, not a modified transition generator. The original U may remain unresolved even when future rank responses are strong. The single selected lag does not estimate a dynamic exponent or support continuum field counting. All new views share their original permutations and covariance.
