# Positive safe-site Gamma does not mean that local loops promote birth

A positive first/completion Gamma among rank-preserving next sites excludes
mixing with an immediate-birth site, but **does not exclude mixing between
different safe mechanisms**. An explicit torus example below has positive
Gamma while adding a contractible loop reduces both near-future birth
responses relative to a tree extension. The role-resolved fork can distinguish
these possibilities without dividing by noisy per-prefix role frequencies.

## Exact role decomposition and pair masks

Fix the ordered prefix Z and let `m(u)=E[X|Z,u]` be the complete future
response. A discrete next-site role C(Z,u) must be determined before the
suffix. Write `pi_c=P(C=c|Z)`, `mu_c=E[m|Z,C=c]` and
`B_c=Cov(m|Z,C=c)`. Then

\[
 B=\operatorname{Cov}(m|Z)
 =\sum_c\pi_c B_c+\operatorname{Cov}_c(\mu_c|Z).
\]

Use the unchanged nested estimator from `84018e19`,
`Dhat=(ab^T+ba^T)/4`, where a,b are the two independent-suffix differences
between iid next labels U,V. Its conditional expectation is
`(m(U)-m(V))(m(U)-m(V))^T/2`. Therefore

\[
 E[1_{C(U)=c,C(V)=c}\widehat D|Z]=\pi_c^2 B_c,
\]
\[
 E[1_{\{C(U),C(V)\}=\{c,d\}}\widehat D|Z]
 =\pi_c\pi_d[B_c+B_d+(\mu_c-\mu_d)(\mu_c-\mu_d)^T],\quad c\ne d.
\]

The mixed mask is unordered; an ordered c,d mask has half its value.
In particular, **safe-safe gives pi_safe squared times within-safe covariance**,
whereas a mixed safe/birth mask contains both within-role matrices as well
as the between-role displacement. It is not a pure gate-between estimator.
The masks sum back to B. Applying the signed projection
`Gamma=(B_AA-B_EE)/4=B_F1,F2` gives the corresponding exact scalar identities.

Safe sites may then be divided into forest attachments and contractible-loop
attachments. Their within-safe Gamma is itself the sum of within-subrole
Gamma and covariance of subrole means. Positive safe-safe Gamma alone does
not identify which of those terms carries it.

## An unbiased common-weight within-role contrast, without pi division

For a quartet j define `M_c(j)=1{both its labels have role c}`. The raw
`M_c(j)*Dhat_j` already estimates `pi_c^2 B_c`, with no normalization.
For a fairer c-versus-d contrast, use two independent quartets at the same
prefix. With M quartets, average

\[
 \widehat T_{c,d}=
 \frac1{M(M-1)}\sum_{j\ne k}M_c(j)M_d(k)
 (\widehat D_j-\widehat D_k).
\]

Its exact target is

\[
 E[\widehat T_{c,d}|Z]=\pi_c^2\pi_d^2(B_c-B_d).
\]

Both roles now have the same nonnegative prefix weight. The Gamma projection
directly tests their **within-role** first/completion coupling contrast.
It needs only the already independent label quartets and two suffixes per
label, not fitted means or a small estimated pi in a denominator. Across
prefixes its target is `E[pi_c^2*pi_d^2*(B_c-B_d)]`, not an unweighted
average of conditional covariances. Role rarity can reduce power; the formula
does not eliminate that fact. Overlapping quartet pairs are not independent
replicates: average within prefix and preserve the original20 batches.

## Why safe loop and forest roles can have different future clocks

Let e be the number of occupied contact edges and c the number of distinct
old components touched. Adding a site changes the ordinary graph cycle rank
by exactly `e-c`. If e=c, the attachment is a forest operation: c=0 creates
an isolated vertex, c=1 extends a component, and c>=2 merges components.
For an **R0-safe** insertion every new cycle has zero ambient winding, so
e-c>0 counts independent contractible cycle additions. In R1, rank-preserving
does not imply zero new winding: additional cycles may run parallel to the
existing essential line. R1 requires the gain mark before calling them
contractible; e-c alone does not establish that interpretation.

The contact theorem `e67d9b90` also explains what rank alone loses. A safe
attachment changes component partition, lift footprint and future contact
addresses while leaving the current homology image unchanged. Merging old
components makes previously independent gauges comparable; extension can
bring a lift footprint to a new seam contact. A local zero-gain cycle need
not make comparable progress toward an ambient birth. None of these facts
fixes the sign of a population response.

### One true torus prefix, two safe sites, different first and second births

Use the 5-by-5 NN torus, period matrix `P=5I`, with exactly these eight
occupied sites (coordinates modulo5):

```
O={(1,0),(2,0),(3,0),(1,1),(0,1),(0,2),(0,3),(0,4)}.
```

This is one induced tree, hence R0. Compare two vacant next labels:

- `u_ext=(4,0)`: e=c=1, a tree extension, still R0.
- `u_loop=(2,1)`: e=2,c=1, closes the local plaquette through
  `(1,0),(2,0),(1,1)`, still R0.

Both children have nine occupied sites and16 remaining labels. With root
`(1,0)`, their common old lift potentials are
`p(k,0)=(k-1,0)`, `p(1,1)=(0,1)`, `p(0,j)=(-1,j)`.
The extension has `p(4,0)=(3,0)`; the loop child has `p(2,1)=(1,1)`.
Apply `P^{-1}[(delta-p)_e-(delta-p)_anchor]` to the contacts of a further
vacant label. The complete nonzero-generator list is:

| Further label | Extension child | Loop child |
|---|---|---|
| (0,0) | span{(-1,0),(0,-1)}, rank2 | span{(0,1)}, rank1 |
| (4,1) | span{(1,0)}, rank1 | zero |
| (1,4) | span{(0,1)}, rank1 | span{(0,1)}, rank1 |
| (4,4) | span{(1,-1)}, rank1 | zero |

Signs of rank1 generators are immaterial to their span. The remaining
multiple-contact candidates are `(2,1),(1,2)` in the extension child and
`(1,2),(3,1)` in the loop child; their address differences are zero. Every
other vacant candidate has at most one contact. This establishes the whole
one-further-label list by the contact theorem, with no exhaustive engine run.

Thus at total occupation10 the conditional first/completion birth responses
are respectively `(4/16,1/16)` after extension and `(2/16,0)` after the loop.
The state `(rank0,occupied_count9)` does not close even the next-step law.
The full canonical curves cannot be identical because these fixed-occupation
Bernstein coefficients already differ; no sign at the production p_ref or
for its integral is inferred from this small example.

There is also an exact sign lesson. Restrict the first next-label choice to
these two safe sites with equal weight, and use those occupation10 responses
as m1,m2. Then

\[
 \operatorname{Cov}(m_1,m_2)=1/512>0,\quad
 \operatorname{Cov}(e-c,m_1)=-1/32,\quad
 \operatorname{Cov}(e-c,m_2)=-1/64.
\]

Each role contains one candidate in this specified comparison, so its
within-role covariance is zero. Positive Gamma is entirely **between safe
roles**, and the loop role has the lower response in both coordinates.
This is a paper counterexample to an interpretation, not a production
population estimate or a universal direction of loop effects.

## Separate intrinsic response coupling from the same-label CRN coupling

For the paired H4-normalized contrast, let delta be delta_cos4 and write
`m_jf(U),m_js(U)` for birth j in the two orientations. At fixed Z,

\[
 \Gamma_{\rm same}=\delta^{-2}
 [\operatorname{Cov}(m_{1f},m_{2f})+
  \operatorname{Cov}(m_{1s},m_{2s})-
  \operatorname{Cov}(m_{1f},m_{2s})-
  \operatorname{Cov}(m_{1s},m_{2f})].
\]

All four terms are available as projections of the same quartet matrix.
The first two are orientation-intrinsic next-label terms. A hypothetical
independent next-label coupling, preserving these same marginal uniform
address distributions, has

\[
 \Gamma_{\rm independent}=\delta^{-2}
 [\operatorname{Cov}(m_{1f},m_{2f})+
  \operatorname{Cov}(m_{1s},m_{2s})].
\]

Their difference isolates the two cross-orientation common-label terms;
it can be read algebraically without a new relabel simulation. A single
fixed permutation relabel does not generally zero those terms. Independence
requires independent labels or averaging an independently randomized uniform
relabel coupling. This reference concerns the next-label Doob layer, not
removal of the common-prefix covariance or of all shared-suffix noise.
For role-conditioned comparisons, keep the marginal allowed label sets
fixed explicitly; changing the role selection at the same time changes the
estimand. A paired Gamma sign alone cannot be assigned to an intrinsic
orientation mechanism without retaining these four components.

Scientific card: the new discriminants are the same-weight role-internal
Gamma contrast and the intrinsic-versus-CRN four-term split. The explicit
R0 tree shows why contractible-cycle count and current rank cannot identify
future first/completion response, and why positive safe Gamma need not mean
loop facilitation. All formulas preserve prefix-measurable roles and the
original full20k denominator. The current production/contact analyses can
read these named quantities with their one shared covariance; no new
sampling, DP, graph enumeration, tests, PR or comment is performed here.
