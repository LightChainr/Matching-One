# Ordinary scalar transport: the N100 inter-peak valley fills in

**The no-Jacobian alternative makes a different, testable prediction.**
N100 A is single-sign but double-peaked. A monotone scalar reparametrization
can move both peaks and the intervening valley; it cannot change their
ordered heights relative to one common amplitude. The existing stream
changes both height ratios, with an exploratory joint Gaussian-reference
`chi2=9.64772/2`, nominal `p=0.00804`.

This is not an extrapolation of the density-warp result. It directly tests
a necessary condition of the ordinary scalar hypothesis on the original
N100 histogram. There is no new sampling.

## 1. Exact identifiability: a whole ordered curve, not pointwise branches

Write `D_j=Y_j(4i)-Y_j(2i)` and
`U_j=Y_j(1/2+i)-Y_j(2i)`. The hypothesis is

\[
U_j(p)=a_jD_j(\phi(p)),\qquad j=A,E,
\]

with independent nonzero constants `a_A,a_E` and one increasing,
endpoint-fixed regular coordinate map phi. There is **no Jacobian**.
Consequently no clock/area ratio is imposed on either amplitude.

For one field, `U'=a phi' D' composed with phi`. If phi'>0, it maps the
complete ordered list of critical points bijectively. Critical types reverse
if a<0, but dividing by a removes that overall sign. Thus if the ordered
critical heights are `h_D=(h_D1,...,h_Dm)` and similarly h_U,

\[
\boxed{h_U=a h_D.}
\]

In particular, for nonzero first-peak height the two nulls are

\[
\frac{h_{U,\mathrm{valley}}}{h_{U,1}}
-\frac{h_{D,\mathrm{valley}}}{h_{D,1}}=0,\qquad
\frac{h_{U,2}}{h_{U,1}}-\frac{h_{D,2}}{h_{D,1}}=0.
\tag{1}
\]

These do not assume linear, polynomial, small, or near-critical phi. They
also hold for arbitrary increasing homeomorphisms, because local extrema
and their values are topological invariants of the ordered scalar curve.

Conversely, for smooth one-dimensional profiles with the same finite
sequence of nondegenerate extrema and strictly monotone intervening
branches, matching the ordered critical heights and endpoint values up to
a determines phi branch by branch:

\[
\phi=(D|_{\text{same ordered branch}})^{-1}\circ(U/a).
\]

At a matched nondegenerate extremum,
`phi'=sqrt(U''/(a D''))`; both neighboring branches join with the same
positive derivative. Endpoint vanishing orders must match for regular
endpoint derivatives. Without that extra condition the statement is about
a monotone homeomorphism. Branch correspondence is fixed once by order:
one cannot independently choose the left or right inverse at every p.

For A and E together, the *oriented parameter curve*
`(D_A(p),D_E(p))` must agree with `(U_A(p)/a_A,U_E(p)/a_E)` after one
reparametrization. For example, A's value at each ordered E zero, divided
by A's first-peak height, is invariant. Those are additional common-map
gates; the present result needs only the simpler A-only necessary condition.
A failure of A alone is sufficient to reject the common A/E scalar model.

This explains a useful contrast with density transport: a single-sign A
always admits a monotone **cumulative density** map after area normalization,
but a single-sign *multipeak scalar* still carries nontrivial height
invariants. Positivity does not imply scalar non-identifiability.

## 2. What changes in the actual N100 profile

Data are the public raw histograms at
[`7b30648`](https://github.com/LightChainr/Matching-One/commit/7b30648be558df0652a7ff22143cc87ed399d042)
in [PR #484](https://github.com/LightChainr/Matching-One/pull/484), with
2,000,000 permutations per shape pair and 200 common batches.
The same counter stream joins all three shapes; they are one dependency
block. P4 normalizations, including their signs, are taken from the source
contract rather than re-estimated.

| A landmark | D position | D height | U position | U height |
|---|---:|---:|---:|---:|
| first ordered peak | 0.3792389 | 0.01141979 | 0.4536278 | -0.00296082 |
| central valley | 0.6035367 | 0.00058429 | 0.6005503 | -0.00105794 |
| second ordered peak | 0.7866146 | 0.00827559 | 0.7289896 | -0.00262169 |

The negative U amplitude is removed automatically by the ratios. A's
amplitude-normalized central valley rises from **0.0512 to 0.3573 of the
first peak**. The two peaks also become more balanced, rather than only
moving horizontally.

| Height ratio | D | U | difference U-D | shared-batch SE |
|---|---:|---:|---:|---:|
| valley / first peak | 0.05116472 | 0.35731524 | +0.30615052 | 0.11437829 |
| second / first peak | 0.72467124 | 0.88546128 | +0.16079004 | 0.07630158 |

The complete 2x2 covariance of the differences is

```text
[[0.01308239, 0.00184953],
 [0.00184953, 0.00582193]].
```

The marginal standardized differences are 2.67665 and 2.10730; their joint
statistic is **9.64772 on 2 df, nominal p=0.0080357**. This is approximately
1%-level exploratory tension with an arbitrary monotone scalar map, not a
claim that the scalar mechanism was rejected in an independent preregistered
experiment. It is notably weaker than some previously tested density
constraints; the mechanisms and their statistical strength are kept separate.

## 3. Error propagation and the two-branch boundary

Each leave-one-common-batch-out calculation removes the same batch from
all three shapes and both orientations, reconstructs both A polynomials,
relocates all three ordered critical points, and recomputes both height
normalizations. The covariance uses `(B-1)/B` times the centered LOO outer
products. The full 18-coordinate covariance and all 200 LOO vectors are
saved, including positions, heights, ratios, and ratio differences.

The critical-height first variation has no first-order displacement term
because A'=0 at an extremum. Nevertheless the calculation actually re-solves
the extrema in every LOO sample. Curvature types remain unchanged and all
LOO extrema remain in the original ordered brackets. Small jackknife bias
estimates for the differences are `(-0.0043125,+0.0021783)`, far below their
SEs; reported point estimates are the direct plug-in values.

An exact integer-Bernstein subdivision certificate proves that **each
empirical mean A polynomial has exactly three interior critical points**.
This precludes quietly choosing some other inverse branch in the fitted
mean curves. It is not a proof of the population's number of extrema.
Sparse tail coefficients are not used to infer physical endpoint orders.

The two landmarks and normalizations were selected after inspecting N100
mean curves, so the nominal Gaussian p-value is not selection-adjusted.
No E-derived p-values, density invariants, or previous clock fits are pooled
with it as independent evidence.

## 4. Direct next prediction, without a scale law

At another scale the ordinary scalar hypothesis still predicts the two
dimensionless differences in (1) to be zero whenever the same ordered
double-peak/valley pattern is present. No transfer of an N100 amplitude,
critical exponent, or polynomial velocity is needed. If the pattern itself
changes, report that change rather than picking another pair of peaks.

The N100 differences above are source-established auxiliary targets for a
new N400 stream, **not** a claim that those nonzero magnitudes must persist
with scale. Their direction suggests an inter-peak filling/peak-balancing
response. Whether that disappears, persists, or reverses is new science.

## Reproduce

```sh
python3 scripts/p267_scalar_clock_transport.py
python3 -m unittest discover -s tests -p test_p267_scalar_clock_transport.py -v
```

The script reads fixed source git blobs, records their SHA256 hashes, and
writes `results/p267-scalar-clock-transport/{score.json,REPORT.md}`. It needs
NumPy and SciPy, uses no server, and does not duplicate the raw production.
The source object must be present locally, e.g. by fetching PR484's branch.

## Scientific card

- **Changed mechanism space:** ordinary scalar reparametrization with free
  A/E amplitudes has its own modest but resolved N100 tension: the A
  inter-peak valley fills in and relative peak heights change. This cannot
  be absorbed by merely moving points along p.
- **Not established:** not an independent confirmation, not a continuum
  field identity, not a state count, and not a density-warp argument.
- **Observer / sector / source / geometry:** ordinary P4[A_top] thermal
  profiles, three N100 modulus contrasts, signed integer-period orientation
  normalization, source 7b30648.
- **Dependency:** all readouts reuse the one N100 stream, seed
  `20260831125401`, offset `267100000000`; no additional Monte Carlo.
- **Next discriminant:** source-frozen ordered ratios at a new scale, with
  both peak branches and the intervening valley retained. Null is zero;
  no untested cross-scale amplitude equality is imposed.
