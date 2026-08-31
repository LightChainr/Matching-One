# A regular local pair interaction: zero Q1 endpoint, nonzero Q activation

**The completed interaction has an exactly invisible Q1 endpoint and a
strictly negative mixed original-U response.** On the fixed original N25
pair, with the canonical singlet coefficient identically one,

\[
\boxed{\partial_\epsilon U|_{Q=1}=0,\qquad
\partial_{\log Q}\partial_\epsilon U|_{Q=1,\epsilon=0}
=-0.04503611397592696.}
\]

The first equality holds on every finite graph where the original
observer and regular common root are defined. The second is the one
completed N25 calculation; it is not a continuum assertion. A concrete
two-site contraction also has a nonadditive activated coupling `13/8`.

## 1. The microscopic alternatives have been separated

| Specified construction | Exact consequence |
|---|---|
| Pure C4-averaged pair tensor K2 | Finite single tangent, but a two-insertion Q1 pole on an actual 8x8 torus occupation |
| Canonical Kreg=K2+K0=average i(I-P1)i† | Entrywise regular; every nonempty finite-network insertion vanishes at Q1 |
| Q activation of this same completed tensor | Nonzero one-mark original-U response and a nonadditive two-mark conditional colour response |

The [8x8 witness](local-pair-two-insertion-geometry.md) consists of four
disjoint occupied paths connecting two vacant sites. Its original rank
is zero. Each single K2 insertion closes to zero, while the double
coefficient is

```
v^18 Q^70 ||K2||^2,
||K2||^2 = Q(Q-3)(3Q^2-9Q+8)/[8(Q-1)(Q-2)].
```

Its Q1 residue is `v^18/2`. Independent activities isolate that
occupation monomial, so no other occupation cancels this coefficient.
This excludes the uncompleted tensor as a coefficientwise-Q1-regular
local interaction family. It does not prove divergence of the summed
homogeneous partition or of U. The old positive
[first-tangent response](local-four-port-transmission-result.md) remains
correct in its stated linear-response scope.

## 2. The completed family carries a genuine two-mark Q interaction

The [kernel algebra](local-pair-two-insertion-algebra.md) gives

```
||Kreg||^2 = (Q-1)(3Q^3-12Q^2+20Q-24)/[8Q(Q-2)].
```

On the same four-path occupation, normalize the marked colour weight
by its unmarked value, and use independent local parameters lambda_x,
lambda_y. The complete first Q derivative of its logarithm is exactly

```
Xi_A = lambda_x + lambda_y + (13/8)lambda_x*lambda_y.
```

This [nonadditivity result](regular-pair-two-site-q-susceptibility.md)
excludes replacing the completed family by independent products of its
single-insertion weights. It concerns a fixed occupation's colour
contraction; 13/8 is not a measured size exponent or a global connected
correlator amplitude.

## 3. One prescribed transmission into the original observer

The [frozen source and full mixed-derivative interface](regular-pair-activation-original-u.md)
come directly from the 15 exterior-component contractions. Because its
relative first coefficient is identically zero at Q1, only its explicit
Q derivative a survives in the mixed response. The original normalizers,
thermal derivatives, induced common-root movement and denominator response
are still all retained.

The exact rational enclosure of `W/A25` is strictly negative:

```
lower = -4818794335234287222311009716950115569 / 10000000000000000000000000000000000000000
upper = -2409397167617143611155504858475047623 / 5000000000000000000000000000000000000000
```

| Full original-U term | Contribution to W |
|---|---:|
| centered direct thermal response | -0.014619566648382092 |
| common-root motion | -0.04875130947142392 |
| source change of thermal slope | -0.03049360208575359 |
| root change of thermal slope | +0.048828364229632645 |
| **sum** | **-0.04503611397592696** |

The mixed root motion is `d_logQ d_epsilon p0=+0.012647482672840315`.
No root-related term is discarded. The exact identity
`W[c0+c1*K]=0` removes any explanation by pure temperature
reparametrization: the source decomposition `a=(1-K/N)+b` gives `W[a]=W[b]`
algebraically, without another score.

The full source coefficients, normalized crossmoments, root data and
rational bounds are in [score.json](../results/p337-regular-pair-activation/score/score.json).

## 4. Scope, execution and the completed decision

Regularity within `K2+c(Q)K0` forces `c(1)=1`, not `c'(1)=0`.
For `c(Q)=1+alpha*(Q-1)+...`, the mixed response changes by
`-alpha*V_old`. No alpha was chosen or fitted. The reported W belongs
to the canonical constant-coefficient completion and cannot be renamed
a completion-independent observable, a Jordan identification, or the
mechanism of the asymptotic Matching-One anomaly.

The contract `25f70f68` and evaluator `5bff5008` preceded the missing
crossmoment collection. Algebra `235f9b0b` and geometry `05d8151f`
were accepted before GO; their public equivalents are `7e46c74c`
and `dcabe53d`. Producer `21d901f9` is public as `3336bf31`;
raw `02fe40f1` is public as `aede8468`. The score ran at `aede8468`;
the sole later pre-score code edit, `60a7f481`, only clarified the
regular-root domain of a metadata sentence, with no numerical change.

Each geometry had one new fixed-source traversal of the same full
`2^25` occupation population: 1.18082 and 1.1344 seconds on this Mac.
The single rational score took 0.244 seconds, importing the old original
root, slope and U/A. [The raw receipt](../results/p337-regular-pair-activation/run.json)
records commands and hashes. No MC, cloud job, root search, old-source
rescore or parameter scan was run. These are new prescribed moments on
the same exact population, not independent stochastic evidence.

The fixed canonical mixed-response zero null is rejected. This finite
question is complete; no c/Q/support/size choice is added to improve it.
Any next field-level claim needs a prediction about occupation-summed
spatial correlations or scaling of a specified completed interaction,
not another same-N25 mixture. Retired P154/P334/F4 production remains
stopped, and no task priority or merge decision is changed here.
