# Regular one-colour local deformations are thermal reparameterizations

**Theorem.** At Q=1, a homogeneous finite site-local tensor deformation
using only the original edge colours, regular under singleton-colour
specialization and keeping the original vacant/occupied summands,
reduces to two scalar site weights. Its normalized occupation law is
iid Bernoulli with an effective density. If those weights are the same
for the two geometries, the original moving-root, slope-normalized U
is exactly independent of the deformation wherever the root and
thermal coordinate are regular.

The tensors need **not** vanish on the all-one colour assignment.
Consequently a regular completion inside this class cannot preserve
the old nonzero direct Q=1 pure-pair U tangent. This is a no-go for
that declared microscopic class, not for all Q=1 theories.

Base: `2ba8863f`. The canonical completion and its finite-network zero
are in [two-insertion algebra](local-pair-two-insertion-algebra.md).
The scalar reduction below is broader than that zero theorem. No new
samples, partition enumeration, fitting or continuum assumptions enter.

## 1. Exact singleton-colour reduction

Fix a finite honest square quotient and its original occupied subset
A. At each site use the same two local tensors

\[
 T_0(\boldsymbol c;Q,v,\epsilon),\qquad
 T_1(\boldsymbol c;Q,v,\epsilon),
 \tag{1}
\]

where the first is assigned to a vacant site and the second to an
occupied site; the activity may already be included in T1. They only
use the original incident edge-colour indices. Assume a finite equality-
diagram expansion with scalar coefficients regular at Q=1, or more
generally an equivalent well-defined singleton-colour specialization.
There is no unresolved colour pole or extra state surviving at Q=1.
The ambient-rank factor is the stipulated `Q^(-r(A)/2)`, and q/E
continue to refer to the original A, not to virtual tensor wiring.

Every edge then has exactly one colour. Set

\[
 w_s(v,\epsilon)=T_s(1,1,1,1;1,v,\epsilon),\qquad s=0,1.
\]

For fixed A the colour contraction is the product of these scalar
values. The rank weight becomes one. Thus

\[
 Z_1(v,\epsilon)=\sum_A w_0^{N-|A|}w_1^{|A|}
 =(w_0+w_1)^N,
 \qquad
 \mathbb P_1(A)=p_{eff}^{|A|}(1-p_{eff})^{N-|A|},
 \quad
 p_{eff}=\frac{w_1}{w_0+w_1}.
 \tag{2}
\]

For a probability statement take real positive w0,w1; near the usual
`w0=1,w1=v>0` this holds for a sufficiently small real deformation.
The partition identity itself is analytic wherever the normalization
is nonzero. Possible common analytic partition prefactors cancel.

This proof uses finiteness to perform specialization before any
infinite-volume limit. It does not replace the Q derivative of a
colour sum by a singleton-colour derivative: those are different
operations, and the new kernel below exploits that difference.

## 2. The original U is unchanged, not merely small

Let M0(p) be the separately normalized, pooled original q expectation
of the fixed geometry pair. Let Y0(p) be its specified angular
difference of E expectations. The same scalar weights at every site
and on both geometries give

\[
 M(v,\epsilon)=M_0(p_{eff}),\qquad
 Y(v,\epsilon)=Y_0(p_{eff}).
 \tag{3}
\]

Follow a simple pooled root `M0(p0)=0` with
`partial_v p_eff != 0`. At its local continuation `v0(epsilon)`,
`p_eff(v0(epsilon),epsilon)=p0`, and

\[
 \boxed{U(\epsilon)
 =A_N\left.\frac{\partial_vY}{\partial_vM}\right|_{v_0(\epsilon)}
 =A_N\frac{Y_0'(p_0)}{M_0'(p_0)}.}
 \tag{4}
\]

The common thermal Jacobian cancels exactly. Therefore all direct
epsilon derivatives of U vanish on this regular root branch, not
only its first derivative. At fixed bare v away from the moving root,
the expectations may change; equation (4) is the original root-tracked
functional, not an assertion of off-root invisibility.

The infinitesimal occupation score makes the same mechanism explicit:

\[
 \partial_\epsilon\log W(A)
 =N\partial_\epsilon\log w_0
   +|A|\partial_\epsilon\log(w_1/w_0).
 \tag{5}
\]

It contains only a constant and total occupation. Normalization removes
the constant and the original U quotient removes the common density
clock. The canonical `1+epsilon Kreg` vacant tensor is a stronger
special case: its all-one insertion is zero, so w0=1 and w1=v and
even the effective density is unchanged.

## 3. The boundary of the no-go

Equation (4) does not automatically cover:

* position-dependent or geometry-dependent scalar weights: these
  generally produce an inhomogeneous law or two different thermal maps;
* a new rank/observable assignment based on inserted diagram joins;
* extra spins, edge states, or direct interactions between neighbouring
  occupation labels that remain nontrivial at Q=1;
* genuinely multi-site occupation factors rather than products of
  the two site-local summands;
* a reassignment of what counts as occupied depending on a neighbouring
  summand or virtual connectivity;
* singular Q coefficients or external projectors, or an infinite-volume
  limit taken before the stipulated finite Q=1 specialization.

These are exclusions from the theorem, not established mechanisms that
recover the previous response. A constant scalar tensor at Q=1 is
allowed and simply changes w0 or w1; making it nonzero does not evade
the theorem. A regular completion retaining the old nonzero direct
tangent must change at least one substantive assumption above.

The bounded finite occupation tilt `exp(-epsilon t(A))` is consistent
with this conclusion. If t is the old component-dependent local-pair
occupation statistic, it is nonlocal as a function of A and is not
an affine function of |A|. It therefore does not factor into common
vacant/occupied scalar weights. This defines a different finite
occupation interaction, not a counterexample to the site-local
edge-colour reduction. If t were affine in |A|, (4) would apply again.

## 4. The new spatial kernel is the Q derivative of log Z

Now use the **canonical** Kreg completion, with distinct site
parameters lambda_x and lambda_y. Let beta_x(A,Q) be its relative
single-insertion colour contraction and beta_xy(A,Q) its relative
double-insertion contraction, including the requisite vacant-site
indicators. These are contractions in the same original occupation
A; beta_xy is not defined as the product beta_x beta_y.

At zero insertion strengths, with expectation in the undeformed
rank-weighted occupation family,

\[
 C_{xy}(Q):=
 \left.\partial_{\lambda_x}\partial_{\lambda_y}\log Z\right|_0
 =\langle\beta_{xy}\rangle_Q
   -\langle\beta_x\rangle_Q\langle\beta_y\rangle_Q.
 \tag{6}
\]

The canonical finite-network theorem gives pointwise, for every A,
`beta_x(A,1)=beta_y(A,1)=beta_xy(A,1)=0`. Define

\[
 a_x(A)=\left.\partial_Q\beta_x(A,Q)\right|_1,
 \qquad
 g_{\pi_{xy}}(A)=\left.\partial_Q\beta_{xy}(A,Q)\right|_1.
\]

Differentiating (6) therefore yields the exact finite-source identity

\[
 \boxed{\left.\partial_Q C_{xy}(Q)\right|_1
 =\mathbb E_{\mathrm{iid},\,p=v/(1+v)}
       [g_{\pi_{xy}}(A)].}
 \tag{7}
\]

The derivative of the background law multiplies the pointwise zero
beta_xy, and the derivative of the one-point product is zero because
both factors vanish at Q=1. The same conclusion holds along the
declared activity path: the extra background/activity derivative again
multiplies the zero endpoint. In particular (7) is **not**
`Cov_iid(a_x,a_y)`. The product of the two first-Q one-point amplitudes
starts at order `(Q-1)^2`, whereas the double colour closure can
already start at order Q-1. This is why the new spatial readout requires
the common eight-port closure rather than recycling two Bell4 scores.

## 5. A cut-factorization support theorem for that kernel

For fixed A, count exterior hypergraph components that contain ports
of **both** x and y. If this count is at most one, then

\[
 \beta_{xy}(A,Q)=\beta_x(A,Q)\beta_y(A,Q),
 \qquad g_{\pi_{xy}}(A)=0.
 \tag{8}
\]

To prove the first identity at integer Q, regard the independent
uniform colour of each exterior component as the random variable.
With no shared component, the two local contractions factor. With one
shared component, condition on its colour c. All other component
colours used by the two sites are independent. The conditional colour
sum at each site is a function of the single c; global colour
permutation invariance makes that function constant in c. Averaging
over c consequently leaves the same factorization. Multiple ports of
one site may belong to this shared component; the argument is
unchanged. Rational continuation proves the identity for the declared
diagram family. The second equality in (8) follows from the two
pointwise zero single-insertion endpoints.

Thus at least two shared physical exterior components are necessary
for this first-Q spatial kernel. They are not sufficient for a
nonzero value or a prescribed sign. The four-path 8x8 construction
already supplies a nonzero example: its canonical double contraction
has derivative 13/8, from the known Kreg norm, without changing any
occupation or recomputing a population. The mechanism in (7)-(8) is
shared-component colour activation, while the direct homogeneous
Q=1 deformation remains only a density clock as proved in (4).
