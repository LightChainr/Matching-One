# The decimation-forced plaquette source reaches the original global U

## Finite mechanism result

The fixed N25 Gaussian pair gives **V_F4=0.194414686461** in bulk source units.
The exact rational enclosure of V_F4/A excludes zero.
Decision: `F4_thermal_alias_and_bare_cluster_endpoint_U_transport_rejected`. This compares the forced decimation operator with
a thermal-only alias in the same original U, not a new fit or a lag1 rescue.

## Source and observer held fixed

All 2^25 configurations of each quotient (5,0) and (4,3) were enumerated: 67,108,864
configurations in total. Sources are Ctot=CB_NN+CW_matching and the number F4 of
fully occupied unit faces, with normalized weights exp(t*S) and no density factor.
q=CB-CW-(K-T_NN+F4), E=q^2 follow digital Alexander on these honest unit-cell tori.
P4 uses the exact direction difference 1152/625. The root is the new pair's pooled
matching root, not an old production calibration; A=25^(13/8)/2.

## Exact coefficients give a visible endpoint correction

| Quantity | Value (numerical evaluation of exact coefficients) |
|---|---:|
| Pooled root | 0.592665539328227 |
| Native U25 | 0.880466156963 |
| V25 for bulk cluster source | -0.401160031467 |
| V25 for bulk F4 source | 0.194414686461 |
| Bare-cluster prediction for N50 endpoint | -1.2373466865 |
| Complete forced-source prediction for N50 endpoint | -0.637689818343 |
| Missing endpoint correction | 0.599656868157 |

The endpoint dictionary is C_parent=C_child+F4. Consequently the difference between
complete and bare predictions is exactly 2^(13/8)*V25_F4, including root and slope
motion. Parent generators are (5,5),(1,7); the complement sign and period rotation
cancel in U. These are theorem-transported endpoint derivatives, not parent
endpoint simulations. A nonzero F4 term makes the configuration-level failure of
bare source closure visible to the specified global observable.

## Calculation includes normalized covariance and the moving root

For each geometry j_q=Cov(q,S), j_E=Cov(E,S). Define Q=mean(q), Y=P4(E),
D=Q_p, r=Y_p/D. The evaluated response is

`V_S/A = jY_p/D - Y_pp*jQ/D^2 - Y_p*jQ_p/D^2 + Y_p*Q_pp*jQ/D^3`.

The four terms, source root shifts and full rational bounds are retained in
`latest.json`. Per-K integer sums already include binomial multiplicity.
All thermal derivatives include the derivatives of covariance centering.
The physical root is bracketed by exact rational bisection; rational interval
arithmetic then encloses the reduced response. Positive irrational area factors
are applied only for numerical presentation and do not change zero exclusion.

## Scope and uncertainty

No sampling errors or confidence levels apply to this exhaustive calculation.
The rational bounds propagate the root interval conditional on the supplied exact
graph counts. The two quotients have different Smith classes (Z5xZ5 and Z25);
the result is a finite-pair mechanism counterexample, not an asymptotic H4 law,
the N65/N85 production family, a continuum field identity or interior saturation
curvature. The prior P154/P334 stop decisions remain unchanged.

## Subsequent exact source completion and next question

The accompanying [closed-source proof](../../../notes/decimation-closed-source-and-global-u.md)
now completes the dictionary: F_parent=T_child-4K_child+2M and
T_parent=4M-4K_child. Thus S_hat=C+F+T-4K+2N is exactly unchanged by this
endpoint decimation. This forced finite source family closes without a fitted
correction or another descriptor. The note contains the proof; this numerical
calculation establishes its first otherwise-missing F4 contribution to the
specified global U.

An interior transmission law for this same closed source remains open. The
endpoint identity alone is not that law, and no V_T value or interior curve is
claimed here. Repeating this F4 calculation or reopening a failed P154/P334
parameterization is not the next target. No new production block is launched.

## Source and reproduction

The dictionary is pinned at execution commit 207436518db46dd13ef0ec91168cb1c99d52eaea,
`notes/p337-checkerboard-decimation-global-u.md`; the topology proof is
56838d5f068f6f0ba7795926dc9343229bdd28ce, `notes/square-checkerboard-endpoint-homology.md`.
The contract and both scripts are pinned by `code_commit` in `latest.json` and
hashes in `run.json`. The stored score predates the subsequent source-completion
note; that explanatory addition does not change its counts or numerical result.
The table is used for exact value lookup rather than an inferred trend.

Run `python scripts/analyze_decimation_plaquette_u.py --output-dir NEW_DIRECTORY`.
Add `--counts-dir results/decimation-plaquette-u` to consume the saved integer
profiles without enumerating configurations again. Existing outputs are not overwritten.
