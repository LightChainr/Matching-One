# P34: Euler and motif control variates

This item supersedes wrapping-only GLS (C03). The C03 negative result is
preserved as provenance: on this matching construction the five matching-odd
wrapping channels are configuration-identical, so frozen GLS weights are
exactly uniform and the variance-reduction ratio is exactly `1.0x` at
`N=65,85`. That GPU gate failed for a structural reason, not lack of samples.
P34 replaces those duplicate channels by Euler and local-motif controls with
analytically known means.

## Exact identity

On the implemented quotients the configuration-level identity is

```text
C_black - C_white = q + V - E + F0
```

with `q` the common wrapping-difference topology variable in `{-1,0,+1}`,
`V` the occupied-site count, `E` the occupied covering NN-edge count (`2N`
edges), and `F0` the occupied elementary-face count (`N` square plaquettes).
Centering by the Bernoulli means `N chi(p)`, `chi(p)=p-2p^2+p^4`, recovers
the specified equal-mean form

```text
D_cluster = q + (V - N p) - (E - 2 N p^2) + (F0 - N p^4).
```

Exhaustive enumeration:

| geometry | N | configs | identity | wrapping channels |
|---|---:|---:|---|---|
| axis L=2 | 4 | 16 | PASS | identical, q in {-1,0,+1} |
| axis L=3 | 9 | 512 | PASS | identical, q in {-1,0,+1} |
| gaussian (2,1) | 5 | 32 | PASS | identical, q in {-1,0,+1} |
| diamond L=2 | 8 | 256 | PASS | identical, q in {-1,0,+1} |

The C++ self-test additionally exhausts Gaussian `N=5` and `N=13` and checks
the same Euler identity against the production union-find. Monte Carlo
production recorded `identity_l1=0` on every replica at `N=65` and `N=85`.

Conditional motif means given occupancy `K=k` match the falling-factorial
formulae

```text
E[E | K=k]  = 2 N (k)_2 / (N)_2
E[F0 | K=k] = N (k)_4 / (N)_4
```

and the analogous expressions for diagonal dimers, straight 3-paths, and
plaquette corners. Centered residuals are exactly zero on every tiny
`(geometry, K)` cell (path-3 motifs are omitted where the three sites are
not distinct, e.g. axis `L=2`).

## Duplicate wrapping channels

GLS over `D_cross, D_both, D_either, D_direction_0, D_direction_1` is
rejected. `FrozenEstimator.fit` raises `DuplicateChannelError` on the five
identical wrapping differences, both on exhaustive tiny tori and on the
production batch files. No ridge pseudo-inverse is applied.

## Monte Carlo protocol

- Bernoulli site occupation at frozen `p_ref=0.592746050790`
- same-N pairs `N=65` `(8,1)/(7,4)` and `N=85` `(9,2)/(7,6)`
- seed `20260834`, replica counters `[0, 1000000)`, 50 equal batches, 8 OpenMP threads
- SplitMix64 counter RNG keyed by `(seed, N, replica, cyclic site)`; the two
  representations of one `N` share `U_j`
- pilot fraction `0.2` (batches 0–9, 200,000 replicas) freezes OLS weights
- evaluation uses the remaining 800,000 replicas
- canonical centering uses analytic Bernoulli means; microcanonical centering
  uses `E[motif | K]` computed per replica in the C++ kernel

Wall time: `2.05 s` at `N=65` and `2.65 s` at `N=85` (about `4.7 s` combined).
No GPU was started.

## Variance reduction versus best single

`q` (wrapping matching difference) is a much lower-variance single estimator
than `D_cluster` (cluster-count matching function). Best single is therefore
`q` at both sizes. Frozen Euler/motif OLS is evaluated on the held-out
replicas. Ratios are `Var(best single) / Var(adjusted)`.

| N | orientation | plan | VR vs q | VR vs D_cluster |
|---:|---|---|---:|---:|
| 65 | (8,1) | Euler canonical | 2.271 | 18.28 |
| 65 | (8,1) | Euler + motifs, microcanonical | 2.345 | 18.87 |
| 65 | (7,4) | Euler + motifs, microcanonical | 2.341 | 18.80 |
| 85 | (9,2) | Euler canonical | 2.104 | 23.74 |
| 85 | (9,2) | Euler + motifs, microcanonical | 2.164 | 24.41 |
| 85 | (7,6) | Euler + motifs, microcanonical | 2.161 | 24.40 |

The identity rewrite `D_cluster = q + Z_V - Z_E + Z_F` is exact but is not
the minimum-variance combination: OLS coefficients differ from `(1,-1,1)`
and beat both `q` and `D_cluster`. Extra 2x2 / short-path motifs add a
further ~3% reduction. Microcanonical centering is essentially tied with
canonical centering at these sizes.

Orientation-difference OLS from 10 pilot *batch means* overfits extra
motifs (`VR < 1`) and is not used for the gate. Replica-level `Δq` would
need a joint first/second gram, which was not accumulated; `ΔV=0` exactly
under the shared cyclic occupancy, as expected.

## GPU 2x gate

The gate is `VR >= 2x` versus best single on multiple sizes, using the
prespecified Euler-plus-motifs microcanonical estimator:

```text
N=65 (8,1): 2.345x  PASS
N=65 (7,4): 2.341x  PASS
N=85 (9,2): 2.164x  PASS
N=85 (7,6): 2.161x  PASS
gate: PASS on two sizes
GPU started: no
```

A later GPU kernel may therefore carry these controls. This CPU item does
not start that work.

## C03 provenance

The wrapping-only GLS result on this branch remains: five matching-odd
channels identical configuration-by-configuration, GLS VR exactly `1.0x`.
See `results/server-20260828/C01/REPORT.md` and the C03 paragraph of
`results/server-20260828/REPORT.md`. P34 does not rerun that wrapping-only
combination except to reject it.
