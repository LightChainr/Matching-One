# A strong cluster-source response does not yet identify directional H4

## Million-sample production resolves the source coupling, not its direction difference

The existing P40 production Gram matrices now give an actual
**absolute-cluster fugacity response**, without generating or replaying a
configuration. At both N65 and N85, each orientation's matching response
is strongly negative for the raw cluster source and strongly positive for
the separately clock/Euler-compensated source. These are two different
specified perturbations, not a sign reversal of one fixed coupling.

The orientation difference remains unresolved: after full controls its
z-scores are **−.189 and1.360**, with nominal two-size p=.38964.
Thus strong coupling to the global matching readout does not automatically
explain the observed directional H4 correction. This completes the
million-sample q/source task; it is not another pilot-readiness result.

## The physical source and the archived measure

Write `q=rank_black−1`, `S=(c_blackNN+c_whiteMatching)/N` and
`J=Cov(q,S)`. This is the lambda derivative under the finite positive
measure proportional to `Bernoulli(p) exp(lambda S)`, at
`p=.59274605079`. White components use the matching graph; inactive
union-find sites are not components. S is normalized by N, so coupling
to the extensive count would multiply J by N.

For each N, two orientations share each site's random occupation:
N65 uses (8,1)/(7,4), and N85 uses (9,2)/(7,6). Each size has one million
configurations in100 aligned batches. All direction projections below are
`(first−second)/Delta cos(4theta)`, using the exact rational geometry.
No power law, free exponent or continuum field is fitted.

The named **clock** span is `(K,K(K−1))`; the **full** span adds occupied
NN edges T and local Euler count chi. Linear projections include an
intercept by centering and are fitted separately for each geometry.
The full compensated source has zero empirical covariance with these
four controls. At the population level the corresponding covariance
projection is the first-order source direction holding their means fixed.
This is not exact conditioning on K or an RG thermal subtraction.

## Strong orientation-wise coupling survives the controls

For N65, the first orientation has raw `J=−.0037431545 ±.0000210394`,
while the full compensated source gives `+.0018788669 ±.0000109686`.
The second gives raw `−.0037311251 ±.0000208806` and full
`+.0018817820 ±.0000091883`.

For N85, the corresponding first-orientation values are
`−.0031392290 ±.0000175497` and `+.0016333975 ±.0000092425`;
the second gives `−.0031349554 ±.0000195292` and
`+.0016147138 ±.0000095831`.

Every full-compensated orientation-wise response exceeds160 standard
errors in this production analysis. The source retains approximately
65.00–65.09% of its variance at N65 and67.32–67.57% at N85 after the
declared full projection. Consequently, the effect is not confined to
that linear clock/Euler span. It does not follow that the residual is a
local canonical energy operator, that its coefficients are universal, or
that these data identify a second continuum field.

## The small H4 difference is not another source-detection problem

The raw directional J is `−.00000882359 ±.00001971007` at N65 and
`−.00000268032 ±.00001619300` at N85. Their nominal zero statistic is
`.227806/2`, p=.89234.

After the clock controls it is `+.00000369369 ±.00001331793` and
`+.00000726340 ±.00001098375`; nominal p=.77328.
After all controls it is `−.00000213821 ±.00001132698` and
`+.00001171787 ±.00000861650`; nominal statistic `1.88505/2`, p=.38964.

These intervals describe the remaining directional resolution. They are
not evidence of exact zero, and the three correlated projections are
not independent model votes. N85 is not the N130 child, so this result
neither confirms nor refutes the previous N130 auxiliary2.47-SE hint
at the same target. The old20k N65 estimate and this million-sample N65
estimate also remain separate retrospective outputs, not a rewritten pilot.

## The stored paired Gram exposes a cancellation, with a coupling boundary

The cross-geometry Gram permits the exact identity

```text
Delta Cov(q,S) = Cov(Delta q, (S_first+S_second)/2)
              + Cov((q_first+q_second)/2, Delta S).
```

It applies to the raw source and to the separately compensated sources.
For the full projection at N65, the normalized terms are
`+.00008130310 ±.00000763811` and
`−.00008344131 ±.00000879968`; their sum is the small directional J.
At N85 they are `−.00005072097 ±.00000646788` and
`+.00006243884 ±.00000633692`.

The two terms are individually resolved and strongly cancel, but their
allocation depends on the archive's chosen **common cyclic-site coupling**.
They must not be promoted to two independent fields, causal mechanisms,
or additional independent H4 detections. The diagonal-response difference
is a marginal property; this allocation uses cross-geometry joint data.
It is an exploratory explanation of this coupling, beyond the declaration's
primary marginal-response question, not a preregistered mechanism score.

## Source reconstruction and uncertainty

The input is commit `291854a518b03eef4293431b89254f0f4429da53`,
`results/local-20260830/P40-production-motif-projection/N{65,85}/mc.motifs.jsonl`.
The metadata and20-variable Gram definitions were read from that commit.
The existing motif `E` is an **edge count**, not E_top; `F0` counts fully
black square faces. The exact definition maps to

```text
K = V
T = E_edges
chi = V-E_edges+F0
K(K-1)/[N(N-1)] = (E_edges-E_mc)/(2N).
```

The implementation transforms the stored sums and within/cross-direction
Gram matrices by this fixed linear map. Source covariances use the pooled
sample denominator `samples−1`. Each delete-one removes the same original
batch from both directions and refits both projections, preserving27
outputs' complete covariance and all100 delete-one vectors. Floating-point
Gram arithmetic and nominal Gaussian errors are not exact rational or
finite-sample confidence certificates. No redundant full covariance is
inverted.

One earlier navigation inference is corrected by the actual generator:
`counter_uniform` includes `splitmix64(N)` in its key. Equal seed and
counter ranges therefore do **not** mean N65 and N85 share site uniforms.
They are distinct N-domain streams, treated as independent in the nominal
two-size statistic under the usual PRNG-domain model. Within a size, both
orientations and all derived views share their original block. The earlier
metadata-only dependence warning is preserved as history, not used as a
restriction on this analysis.

## Next: distinguish global source coupling from anisotropic coupling

The concrete next physical question is whether this **named compensated
source** couples to E_top and to the norm-4 directional residual along
the actual child/thermal coordinates. Another q-only source-detection run
does not answer it: q coupling is already strong, whereas its directional
difference is much smaller at these two parent geometries.

The archive does not contain `q²*S` or all `q²*control` moments, so it cannot
supply that E_top response by squaring an average or relabeling `E_edges`.
Nor does it contain the mixed third moment needed for the thermal
derivative of J. The existing Phase-E replay did obtain E_top/source at
N65/N130, but at only20k. A next acquisition or old-counter reobservation
should target these explicitly missing mixed readouts with the named source,
not another generic algebra compiler or a repeat of this million-sample q
analysis. First check for a relevant newer result before commissioning it;
this is coordination, not a lock or approval gate.

The cancellation allocation suggests a separate idea: compare physically
defined couplings/lifts if their spatial meaning is the question. Merely
changing cyclic labels can change the two summands without changing the
diagonal H4 response, and is not by itself a new physical experiment.

## Reproduce from the existing source commit

```bash
git fetch origin analysis/p40-production-motif-projection-20260830
python3 scripts/analyze_p40_absolute_cluster.py --output-dir results/p40-absolute-cluster-reproduction
```

The scorer reads immutable Git objects and checks the source hashes in
the [declaration](../../analysis/p40_absolute_cluster_reanalysis.json),
committed at `ae3d55af` before this aggregate. It does not execute the old
Monte Carlo engine or copy a weaker substitute for missing fields.
The [JSON](latest.json) stores definitions, input hashes, coefficients,
paired allocations, full covariance, delete-one vectors and environment.
It explicitly marks the missing E_top and N130 queries `not_scoreable`.

The complete analysis took approximately **.11 seconds on local ARM64
Python**, with zero new samples, zero configuration replays and no scientific
test suite. GPU capacity would not help this small saved-matrix calculation.
The report keeps the technical definitions, evidence and interpretation
separate; two parent sizes are presented as a numerical readout rather than
a misleading scaling plot.
