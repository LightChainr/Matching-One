# Block count does not determine the stationary response

## Technical summary

The fixed block-fugacity prediction fails for the existing width8 P398
join/detach process. Its stationary score has **34.27% relative weighted RMS
error**. More importantly, **10.99% of the true score variance lies within
equal-block-count classes**, so changing the coefficient or allowing an
arbitrary function of block count cannot recover the complete derivative.
Two explicit four-block states even have opposite stationary response signs.

The error matters to the original propagation readout. Substituting the
candidate stationary derivative while retaining the exact saved generator
response predicts a positive +− zero-frequency derivative, **+.01812**;
the actual derivative is **−.05023**. A correct symmetry and a tolerable static
amplitude are not sufficient to determine the integrated response.

The calculation took **.218227 seconds**, using the already saved π, π′ and
186-dimensional character-sector matrices. There was no stationary solve,
new η point, width expansion, finite-difference sweep, Fréchet rerun,
Monte Carlo, server operation or research test suite. Full statewise scores
and readouts are in [latest.json](latest.json).

## One fixed candidate, no fitted coefficient

The model and clock are those of the
[completed η=0 calculation](https://github.com/LightChainr/Matching-One/blob/fb01c44aa45e4f8d37d52144e2ad7c4adfe6ce40/experiments/p398-linear-response-20260831/README.md):
1430 noncrossing partitions, join rate `1+η`, detach rate `1−η`, and sixteen
attempts per time unit. Let B(x) be the number of partition blocks. The candidate
was fixed before this calculation:

```text
pi_eta(x) ∝ pi_0(x) [(1−eta)/(1+eta)]^B(x),
g_candidate(x) = (d/deta) log pi_eta(x)|0 = −2[B(x)−<B>_0].
```

Here `<B>_0=4.5`. The true score is `g=pi_prime/pi_0`, read from the saved
Poisson derivative. Both scores have zero mean. The candidate is correctly
Kreweras-odd, to about1.1e−14; its failure is therefore more specific than a
parity mismatch.

The true stationary tangent equation residual is1.03e−15. Substituting the
candidate `pi_prime=pi_0*g_candidate` yields a maximum residual .00928793.
These are deterministic finite-state discrepancies, not Monte Carlo test
statistics or uncertainty estimates.

## Equal block count hides opposite responses

In restricted-growth-string notation, the following two states both have B=4:

```text
(0,1,1,2,2,3,3,0): pi0=.000322942880, true score=−1.780435265
(0,1,0,0,0,2,0,3): pi0=.000968828639, true score=+3.113768598
```

The block-count candidate assigns **+1** to each. The true score span is
4.894203863. Both states have positive stationary probability. Their block-size
profiles differ, so this witness establishes a need beyond block count; it does
not yet establish a need beyond the entire block-size profile.

The weighted score variance is12.87669987, while the fixed-candidate squared
error is1.512513706. For a stronger sufficiency diagnostic, form
`g_B=E_pi0[g|B]`. The residual variance is1.414804541, or10.9873% of the total.
This is the irreducible error of all block-count-only first-order scores in the
same L² measure. `g_B` uses the true derivative and is only a diagnostic; it is
not an independent prediction or a fitted replacement physical model.

## The missing measure structure changes a response sign

Keep the original two source rays, metric and attempt clock. The candidate
static cross-covariance derivative is `.8839934615`, versus the actual
`.7821171947`. This moderate static-amplitude error becomes much more important
after the stationary and dynamic terms combine.

For the −+ integrated normalized response, the candidate predicts
`−.0077566680` versus the true `−.0530034146`. For +− it predicts
`+.0181196545` versus the true `−.0502337767`: the sign is wrong.
At the old lag t=.5, its +− derivative is already `+.0026756246`, while the
actual derivative is still `−.0680760618`; at t=1 it predicts `+.0233822634`
versus `−.0024588391`. No new zero-crossing search was performed.

Only the stationary metric derivative was replaced. The true saved generator
contribution was held fixed, so this comparison isolates the proposed measure
law. It is not a newly constructed autonomous stochastic model with a different
generator. It also explains why the successful16-column dynamical compression,
which imports the full π and π′, did not already predict the measure mechanism.

## Definitions and reproducibility

In the original character-i basis Q, set
`B_metric=Q* diag(pi0) Q`, `B_metric_prime=Q* diag(pi_prime) Q` and
`H=−Q*F0^T Q`. The source matrix Z, H and baseline metric are read from the
frozen archive. The candidate replaces only `pi_prime` by `pi0*g_candidate`.

The same convention `C(t)=conj[Z* B_metric exp(−Ht) Z]` is used. For the
zero-frequency response `R=C0^−1 conj[Z* B_metric H^−1 Z]`, the stationary part is

```text
Rprime_pi = C0^−1 [conj(Z* B_metric_prime H^−1 Z) − C0prime R].
```

The saved true generator term is then added unchanged. Old-lag stationary
contributions use the same formula with `exp(−Ht)` and the saved baseline U(t).
This avoids another Fréchet or Poisson calculation.

The [contract](../../analysis/p398_block_count_measure_contract.json) and
[script](../../scripts/analyze_p398_block_count_measure.py) were committed at
`7afcbae0260603fb8769ea1b754bae54ebe42234` before execution. Every loaded
`fb01c44a:commit/path`, code hash and local Python version is recorded in JSON.
The pinned geometry constructors are reused without running their analysis
entrypoints. If the source objects are absent in a fresh checkout, fetch
`research/results-first-synthesis-20260831` before running:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python scripts/analyze_p398_block_count_measure.py
```

The saved result is preserved; reproduce in a separate checkout/output copy.
The evidence is a deterministic float64 finite-state computation, not a rational
certificate, statistical replication, continuous-field count or square-site
matching identification.

## Next discriminating question

Separate **block-size composition** from **arrangement within the same block-size
profile** before proposing another one-variable measure law. The present witness
does not distinguish them. A statewise conditional-score comparison on the
existing finite archive can do so without new widths or η values.

This is a finite-model explanation line, not a reason to keep expanding that
model indefinitely. Connecting it to Matching-One still requires an explicit
microscopic source/readout or transfer-operator map. In parallel, P334's new
joint-birth moments call for conditioning on the prefix to distinguish changes
in its internal continuation law from changes among prefix-specific clock shifts.
