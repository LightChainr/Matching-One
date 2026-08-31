# From a resolved even tangent to the original norm-4 source response

The [million-sample result](../results/p40-even-given-odd/REPORT.md) resolves
an absolute-cluster response outside the pure-q tangent of the finite
topological marginal. Its fixed-p quantity C is not the old norm-4 U.
The distinction determines the next useful observation; it does not add
a new approval step or postpone use of the completed result.

## Preserve the actual frozen object

At freeze `2236d36c80c8a466d9317c929bc33e92a7ca9d33`,
[`predictions/norm4_two_generator_transfer_20260829.yaml`](../predictions/norm4_two_generator_transfer_20260829.yaml)
defines `U=N^(13/8) P4_legacy[S_legacy']/Mbar'`.
Since `S_legacy=1−E_top/2` and its direction convention is second−first,
the standard first−second projector gives

```text
U_N = (N^(13/8)/2) B_N/D_N,
B_N = P4_std[partial_p <E_top>],    D_N = partial_p <q>_bar,
all evaluated at the pooled matching root p0(N).
```

The two chains are65→130→260 and85→170→340. The actual residuals in
[`score_norm4_production.py`](../scripts/score_norm4_production.py) are
`U_N−3U_2N+2U_4N` for q2 and `U_N−2U_2N+U_4N` for Jordan.
N130/N170 are intermediate norm-2 steps, not the norm-4 endpoints.
The completed old model scores are not repeated by introducing a new source.

## A common raw source gives an explicit derivative, including root motion

Use the same positive source family
`P_(p,lambda) proportional to Bernoulli(p) exp(lambda S)` with
`S=(C_blackNN+C_whiteMatching)/N`. Keep this definition fixed as p varies.
Let `Jq_g=Cov(q_g,S_g)` and `JE_g=Cov(E_g,S_g)` on geometry g, and use
a bar for the mean over the two orientations. Implicit differentiation gives

```text
p0_dot = −Jq_bar/D,
U_dot = (N^(13/8)/2) [
    (P4[partial_p JE] + P4[partial_p² <E>] p0_dot)/D
    − B (partial_p Jq_bar + partial_p² <q>_bar p0_dot)/D²
].
```

These derivatives are evaluated at the baseline pooled root, with the
source-independent orientation projector. A future observed source
response of the original residual is `L4[U_dot]`, not `L4[C]`.
This formula does not assert that the unperturbed q2/Jordan law remains
valid under arbitrary source deformations; a mechanism must make that
extra prediction before its source response is used to distinguish it.

For an observable O and p-independent source S,

```text
partial_p Cov(O,S) = kappa(O,S,K) / [p(1−p)].
```

Hence the directly missing source products are the q*S*K and E*S*K
moments at the relevant root/thermal coordinate, together with lower
moments already present in this source's finite-p analysis. Existing
threshold histograms can supply baseline root and p derivatives, but do
not reconstruct absolute cluster counts conditional on rank and K.
The new P40 supplement supplies E*S and E*controls at one p only; it
does not supply these triple products or exact transport to a different p.

## Use the finite result now, without assigning it another observable's law

The observed quantity
`C=JE−Cov(E,q)Jq/Var(q)` is the E response in an explicitly q-compensated
measure at one p. It is useful evidence that a named source has an even
topological tangent. It cannot inherit U's13/8 normalization, its q2/Jordan
coefficients, or the separate root-character1/16 rule without a physical
transport derivation.

The original N65/N13020k absolute-cluster CSVs already contain enough
moments for their own fixed-p C and paired uncertainty; no replay is
needed for that supplementary comparison. They do not provide N260 or
the original root-normalized U_dot. A useful next work block can therefore
measure a concrete raw-source U_dot on one available lineage component,
while theory predicts how that same source perturbs the surviving
mechanisms. It need not first build another generic projection framework.

For geometry-adapted residual sources, state which counterterms are fixed
and which depend on p/lambda; additional coefficient derivatives otherwise
enter. In particular, compensating S against Z and then fixing q does not
simultaneously hold every Z mean fixed. Raw S is the unambiguous starting
point for the cross-geometry physical comparison.
