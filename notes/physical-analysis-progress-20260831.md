# Physical response progress after context recovery

This continuation consumed existing results and produced two new scientific
outputs. It did not launch a tool-building or repeated-testing round.

## Norm-4: the missing sector × local-interaction moment is now measured

[`c0880c2`](https://github.com/LightChainr/Matching-One/commit/c0880c297b40699563e8be537e777ac8cd4084c8)
reobserves exactly the old 20k Phase-E replicas at each N65/N130, preserving
their paired directions. It adds the third mixed moment absent from P40's
second-moment archive. [Report](../results/p154-fixed-k-interaction/REPORT.md),
[full paired covariance](../results/p154-fixed-k-interaction/latest.json).

The fixed-K local-edge response to E_top remains unresolved: z=.152/.891,
joint chi-square=.8173/2, nominal p=.66455. Its raw local interaction and
occupation-only second-thermal-score parts are also unresolved. This does
not establish that the effect is zero or entirely thermal. It does complete
one previously unmeasured microscopic-source question on real configurations.

The [source/emission distinction](p154-local-source-versus-connectivity-emission.md)
is useful for choosing the next physical question: disjoint finite
occupation-cylinder functions have zero connected covariance in the
unconditioned product measure, while a local source can still couple to a
global topological event. A connectivity emission or a named source derivative
is a different object. Neither this observation nor the mixed-response result
identifies the original norm-4 secondary field.

## P398: exact eight-mode propagation and useful two-mode compression coexist

The earlier `e38fe76` already completed the positive width-four A/L two-point
matrix; `b35e100` completed its positive anisotropic family. The separately
read `dbd4081` adds continuous-distance amplitude/metric fingerprints for that
same width-four family. Those completed calculations are not reassigned.

[`8f7a587`](https://github.com/LightChainr/Matching-One/commit/8f7a5875157265e32e9db08c5f3991a9b9ddb86e)
keeps the adjacent-pair and singleton emissions and extends the first cyclic
character to width five. In the actual positive 42-state bond-cylinder
measure, the positive-separation two-point Hankel rank is exactly eight;
the inherited two readouts do not close in two dimensions. The rank is
observable in C(1)…C(7), not merely the dimension of an ambient algebra.
[Exact report](../results/p398-physical-two-point/REPORT.md).

Consuming that exact matrix, rather than fitting eight free exponentials,
gives a more useful conclusion. The two largest propagation eigenvalues are
.109393084868 and .0389423547112. Their parameter-fixed residues approximate
the full two-point matrix to .8996% at d=1, .1523% at d=2 and .004050% at d=4.
Yet the two-readout normalized propagator has a 7.8838% semigroup defect at
lags 1/2. [Mode visibility](../results/p398-physical-two-point/mode-visibility.md).

Therefore an effective slow-mode description can be excellent without being
the exact microscopic state. These are finite-width bond results, not eight
CFT fields, not a universal two-field model and not a failure of the underlying
Markov chain. In particular, a small whole-matrix approximation error does
not guarantee a small error for a cancellation-sensitive projected channel.

## Next work should consume these results

The first attention remains identification of the original norm-4 physical
direction, with this measured local-edge row and its explicit limitations
available as inputs. The P398 line now starts after both the width-four
physical calculation and the width-five extension: compare propagation and
emission-dependent slow-mode errors at declared separations, rather than
equating a finite exact rank with the required number of physical fields or
rebuilding a two-point interface. The continuous width-four invariants have
their own model boundary; they are not assumed to hold at width five.

The [current attention board](../docs/NEXT-TARGETS.md) remains the sole ranked
entry. These are scientific suggestions, not sequential permissions. All
existing Issues stay open/unlocked as before, and Draft #267 stays unmerged.

## Execution scope

- P154: deterministic replay of 40k existing replica IDs, two directions
  each; original counts and rank sums match all 400 batch/direction rows.
  New samples: zero. N65/N130 replay times .162/.256 seconds on this Mac.
- P398: one exact width-five point, .169 seconds; then one numerical
  decomposition of its saved matrix, .0053 seconds. No exact-engine rerun.
- Exact rational/cyclotomic statements are separate from complex128 mode
  visibility and estimated-covariance statistical statements.
- No server/tunnel use, new Monte Carlo production, full test suite,
  Issue closure/lock/rename, PR merge, rebase or force-push.
- Code, data, definitions and source hashes are retained with each result.
  Only delivery metadata/links and whitespace receive a closing check.
