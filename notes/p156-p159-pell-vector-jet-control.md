# Pell/Eisenstein full-homology vector-jet square-bond control

Status: **complete; square-bond vector-control gate passed** for issues #156
and #159.  This run replaces the scalar-only question by a probability-vector
value and its first two thermal derivatives.  It does not run N418 or N780.

## Frozen contract

The control uses square-bond percolation at the two first useful Pell
approximants on opposite sides of the Eisenstein modulus:

| design | period matrix | Pell side | N | tau |
|---|---|---:|---:|---:|
| `pell_Dminus2_N30` | `[[6,3],[0,5]]` | -2 | 30 | `1/2+5i/6` |
| `pell_Dplus1_N56` | `[[8,4],[0,7]]` | +1 | 56 | `1/2+7i/8` |

For each design, 240,000 common-random-number replicas were evaluated at
`p=(0.49,0.50,0.51)` in 48 batches with seed `20260901`.  The vector retains
rank zero, rank two, and every observed primitive unoriented winding line.
The central value, symmetric first derivative and symmetric second derivative
are all accompanied by batch covariance matrices.  The N=4 fundamental Pell
quotient is evaluated exactly over all 256 bond configurations and satisfies
the probability/derivative sum rules `(1,0,0)`.

## Result

The three shortest primitive lines were transformed into the real
nontrivial-C3 coordinate `C`, the reflection-odd coordinate `Q`, and the scalar
coordinate `S`, after subtracting the continuum baseline from the central
value:

| design | `C(p_c)` | z | `Q(p_c)` | z | `S(p_c)` | z |
|---|---:|---:|---:|---:|---:|---:|
| N30 / D=-2 | 0.0079343 +/- 0.0007343 | 10.81 | 0.0004078 +/- 0.0008689 | 0.47 | 0.0032440 +/- 0.0010800 | 3.00 |
| N56 / D=+1 | 0.0043203 +/- 0.0007528 | 5.74 | 0.0008841 +/- 0.0007538 | 1.17 | 0.0030577 +/- 0.0009259 | 3.30 |

Both Pell sides therefore recover the same declared organization: a resolved
positive nontrivial C3 component, a reflection-null `Q`, and a smaller but
resolved scalar residual.  The full observed support is identical at both
sizes:

```text
rank0
rank1:(0,1), (1,-2), (1,-1), (1,0), (1,1), (2,-1)
rank2
```

All three probability vectors conserve mass and there are no incompatible
rank-one/invariant failures.

The thermal jet adds a mechanism split that the old scalar statistic could
not provide.  The nontrivial-C3 first derivative is unresolved at this pilot
precision (`z=-0.43` at N30 and `z=1.59` at N56), whereas the scalar second
derivative is negative on both sides and resolved (`z=-3.94`, `-5.29`).  The
nontrivial-C3 second derivative is not resolved (`z=1.22`, `0.40`).  Thus the
critical nontrivial vector offset and the leading thermal curvature occupy
different representation sectors in this control.

## Decision

The square-bond positive-control gate passes.  Any subsequent Pell/Eisenstein
site experiment should freeze the complete homology probability vector and
its thermal jet; returning to a scalar H4 vote would discard resolved
information.  This small run establishes the vector language and covariance
contract without starting the N418/N780 production.

## Execution provenance

The requested Huawei route was attempted first.  NePnUn and HZsCM6 both
reached `Running` but rejected their saved SSH keys with `Permission denied
(publickey)`.  No key reset was authorized or attempted, and each owned tunnel
was stopped; NePnUn returned to `Ready` and HZsCM6 was sent to `Stopping` for
return to `Ready`.  To finish the bounded analysis, the identical frozen
contract ran locally on the 10-core Apple Silicon host with 8 workers.  The
internal measured wall time was 13.899 seconds (shell wall time 14.05 seconds).

## Reproduction

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 \
python3 scripts/pell_homology_vector_jet_control.py \
  --samples 240000 --batches 48 --workers 8 \
  --seed 20260901 --p0 0.5 --h 0.01 --dps 70 \
  --source-commit a9b9bf3b \
  --output-prefix \
  results/local-20260901/P156-P159-pell-vector-jet-control/result
```

Machine-readable results, batch sufficient statistics, runtime output and
checksums are stored in
`results/local-20260901/P156-P159-pell-vector-jet-control/`.

