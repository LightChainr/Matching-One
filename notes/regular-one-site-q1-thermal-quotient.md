# Every regular one-site Q1 colour interaction is thermal in original U

**Theorem.** A homogeneous interaction on one original binary-occupation
site, expressed by four-port equality diagrams with coefficients regular
at Q1, reduces at Q1 to two scalar weights: one vacant and one occupied.
After normalization it is a common Bernoulli thermal reparameterization.
Consequently its complete moving-root/slope U is independent of every
interaction coupling on a regular root branch.

This includes perturbing **both** the vacant and occupied tensors. It
does not require an unequal-pair factor, zero one-insertion closures,
C4 averaging, or a particular singlet counterterm. In this class there
is no completion that retains the nonzero pure-K2 first derivative
`V_old=+.0018155512845251097`. The search for such an entry-regular,
homogeneous, one-original-site completion is therefore finished.

This is an exact finite-network proof, not a new numerical result. It
does not exclude the bounded occupation reweighting, a specified singular
confluent field, multi-site interactions or a Q-activated response.

## 1. Definition of the class

Fix a finite original graph with N binary occupation sites. For a square
lattice site x let its four incident colour ports have indices
`a=(a_N,a_E,a_S,a_W)`. Write eta_x=0 or 1 for its **original** vacancy
or occupation. Allow a vector of interaction couplings lambda and two
otherwise arbitrary tensors

\[
 T_\eta(Q,p,\lambda;a)
   =\sum_{\pi\in\mathrm{Part}(4)}
       c_{\eta,\pi}(Q,p,\lambda)D_\pi(a),\qquad \eta=0,1.   \tag{1}
\]

Here D_pi is the ordinary equality tensor for the indicated port
partition: colours in a block are equal, and distinct blocks are not
forced different. Coefficients may depend nonlinearly on p and lambda,
but admit a finite diagram expansion regular at Q=1 in a neighbourhood
of the coupling point. This precise regularity assumption is stronger
than finding a finite value after one selected exterior contraction.

Use the same T_0,T_1 at every site and in both geometries of a direction
pair. A coefficient cannot depend on neighbouring occupations, the
exterior connectivity partition, K, rank, site position or geometry.
Ordered-port anisotropy is allowed: coefficients of different diagrams
need not be equal or C4 symmetric. Such anisotropy does not survive as
more than a scalar at one colour.

Join colour ports by the graph's ordinary equality wires and close the
finite network. A rank factor such as `Q^(-r_G(A)/2)` is allowed with
its regular branch near Q1; it equals one there and retains the rank of
the stipulated original occupation A. The observable q=r_G(A)-1 and
E=q² also retain that original occupation/rank. Virtual diagram joins
do not alter q/E. No further occupation-dependent factor at Q1 or
singular external projector is included.

The original family is contained in (1), for example
`T_0=1`, `T_1=v(p) D_all4`, `v(p)=p/(1-p)`. More general normalizations
such as `(1-p),p D_all4` give the same normalized law. Near this original
family the two scalar weights below are positive. Algebraic normalized
identities only require a nonzero partition; the Bernoulli probability
interpretation uses positivity.

## 2. One-colour reduction, without a special kernel identity

Define the regular scalar evaluations

\[
 a_\eta(p,\lambda)=\sum_\pi c_{\eta,\pi}(1,p,\lambda).
                                                               \tag{2}
\]

Every equality diagram with l free closed colour components evaluates to
Q^l. Expanding the finite network in (1) produces a finite sum of such
monomials with regular coefficients. Thus continuation to Q1 commutes
with this finite sum and equals evaluation on the one-colour set. Every
delta then equals one, so T_eta equals a_eta, irrespective of its port
partition or the other occupations.

It follows configuration by configuration that the Q1 weight is

\[
 w_G(A;1,p,\lambda)
     =a_0(p,\lambda)^{N-K(A)}a_1(p,\lambda)^{K(A)}.          \tag{3}
\]

The unchanged rank factor equals one. Spectator equality loops also
equal one, so no exterior-connectivity or homology multiplier remains.
An occupation-independent partition prefactor would cancel below.
This proof does not assert that every insertion is zero: either scalar
in (2) may change substantially with lambda.

Summing over all original occupations gives

\[
 Z_G(1,p,\lambda)=(a_0+a_1)^N,
 \qquad \widetilde p(p,\lambda)=\frac{a_1}{a_0+a_1},
                                                               \tag{4}
\]

and hence

\[
 \mathbb P_{G,1,p,\lambda}(A)
    =(1-\widetilde p)^{N-K(A)}\widetilde p^{K(A)}.             \tag{5}
\]

The effective occupation parameter is the same for both geometries.
The graph still determines the original q/E of a configuration, but its
probability is precisely the ordinary independent Bernoulli probability
at the common parameter p-tilde. This holds at finite interaction
strength wherever (4) is defined, not just at first order.

## 3. The complete original-U quotient removes that common parameter

For a same-N direction pair g=a,b, let `m_g^B(s)`, `e_g^B(s)` be its
ordinary Bernoulli expectations of q,E at parameter s. Define

\[
 M_B(s)=\frac{m_a^B(s)+m_b^B(s)}2,\quad
 Y_B(s)=\frac{e_a^B(s)-e_b^B(s)}{\Delta_4},\quad
 A_N=\frac{N^{13/8}}2,
                                                               \tag{6}
\]

with the fixed nonzero exact direction factor Delta4. Equation (5)
implies for the perturbed Q1 family

\[
 M(p,\lambda)=M_B(\widetilde p(p,\lambda)),\qquad
 Y(p,\lambda)=Y_B(\widetilde p(p,\lambda)).                   \tag{7}
\]

Choose a simple pooled Bernoulli root s_star, `M_B(s_star)=0` and
`M_B'(s_star) != 0`. Follow the corresponding root p_star(lambda) with
`p-tilde(p_star(lambda),lambda)=s_star`, on a branch where
`partial_p p-tilde != 0`. The implicit-function theorem supplies this
local branch near the original family. These assumptions are exactly
what is needed for the old thermal-slope denominator to remain regular.
By the chain rule,

\[
 \boxed{U_N(1,\lambda)
 =A_N\frac{Y_p}{M_p}\bigg|_{p=p_\star(\lambda)}
 =A_N\frac{Y_B'(s_\star)}{M_B'(s_\star)},}                   \tag{8}
\]

which is independent of lambda. All coupling derivatives of this U
vanish, including mixed derivatives among any number of regular
interaction coordinates. The result does not depend on the value of
the area exponent: any coupling-independent prefactor in (8) would
have the same invariance. The displayed 13/8 is the actual repository
definition, not a scaling assumption in the proof.

At first order, the cancellation can be seen directly. Write
`c=(partial_lambda p-tilde)/(partial_p p-tilde)` on a regular branch.
Then `jM=c M_p`, `jY=c Y_p`, so the full source response is

\[
 \frac{\partial_\lambda U_N}{A_N}
 =\frac1{M_p}\partial_p
       \left[jY-\frac{Y_p}{M_p}jM\right]_{M=0}=0.            \tag{9}
\]

Expanding (9) is the existing four-term direct/root/slope formula.
The root itself can move: `partial_lambda p_star=-c`. For example,
`T_0=1+lambda`, `T_1=v(p) D_all4` is entry-regular and has nonzero
thermal motion; it still leaves (8) unchanged. Thus the theorem is
strictly broader than proving a particular Kreg insertion vanishes.
Site-average couplings epsilon/N and extensive couplings are both
covered; multiplying an exact zero by N cannot restore a U response.

## 4. Consequence for the completed pair calculations

Execution's canonical `Kreg=K2+K0` has the additional special property
that its one-colour value is zero. Its vacant and occupied scalar weights
are therefore unchanged, giving not only (8) but an unchanged Q1 root
and unchanged fixed-p expectations. This is the stronger special case
already proved in the
[entry-regular completion algebra](https://github.com/LightChainr/Matching-One/blob/7e46c74ce149d5a0a06d1085eb36eebb1bbe6bdb/notes/local-pair-two-insertion-algebra.md).

The present theorem also covers completions whose one-colour value is
nonzero, arbitrary additional regular four-port diagrams, modifications
of the occupied tensor, and nonlinear dependence on homogeneous
couplings. They can only change the two scalar weights in (3). Therefore
none can preserve the old pure-K2 nonzero first U derivative while
remaining in this entry-regular one-original-site class. This conclusion
is not restricted to choosing a coefficient c in `K2+c(Q)K0`.

There is no contradiction with the
[completed canonical mixed response](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/regular-pair-interaction-result.md),
`partial_logQ partial_epsilon U=-.04503611397592696`. Evaluation at one
colour proves a value at Q1; differentiating Q also differentiates
colour-component counts. A function can be constant in lambda at Q1
and have a lambda-dependent first Q derivative. In particular the
completed two-insertion Q susceptibility 13/8 remains compatible with
the theorem; it is not a size exponent.

The already derived finite-counterterm relation
`W_alpha=W_canonical-alpha V_old` concerns that **mixed Q derivative**.
It cannot create a nonzero direct epsilon derivative at Q1. Here W_alpha
retains execution's mixed-response notation; it is not the separate
size-note convention `W_N=N V_av(N)` for the old bounded occupation
tangent. No alpha is chosen or fitted in this theorem.

## 5. Exact scope of the no-go and the routes it leaves open

- **Singular/confluent coefficients.** A tensor with meromorphic diagram
  coefficients whose poles cancel only after selected closures is outside
  (1). That is precisely how the old pure-K2 one-insertion tangent escapes
  the theorem. The theorem does not exclude every conceivable singular
  completion that is regular after a specified set of physical closures.
- **Bounded occupation reweighting.** `exp(epsilon S_av(A))`, with the old
  connectivity-dependent t_x, is a well-defined Q1 law and retains the
  measured linear V. Its weight depends on exterior occupation connectivity,
  so it is not a product of two site-state-only scalars in (3). It cannot
  be renamed an entry-regular one-site colour tensor without another
  representation that satisfies the theorem's assumptions.
- **Multi-site or finite-patch vertices.** A patch with m original sites
  has up to 2^m original occupancy states. Its one-colour weights need not
  factor into a_0^(m-K) a_1^K. For example, an explicit factor depending
  on neighbouring occupations gives an interacting lattice gas, not the
  common Bernoulli reparameterization. Blocking or introducing such a
  factor changes the one-original-site hypothesis; no zero-U conclusion
  for that enlarged class is proved here.
- **Inhomogeneity or geometry-specific couplings.** Site-dependent scalar
  weights give different p-tilde_x, while geometry-specific homogeneous
  weights give different p-tilde_g. Neither is the common scalar map in
  (7). Using those freedoms to fit a direction response is a different
  source contract, not a completion within this no-go class.
- **Changed observers or nonregular global factors.** Recomputing rank
  after virtual joins, changing the q/E observable, inserting singular
  projectors/seams, or retaining an extra occupation-dependent weight at
  Q1 invalidates the assumptions. A zero thermal derivative, vanishing
  normalization, or switch to a different pooled-root branch also lies
  outside the regular domain of (8).
- **Q derivatives and order of limits.** Neither Q activation nor a limit
  taken before finite-network Q1 evaluation is ruled out. The proof is
  finite and contains no interchange with an infinite-volume limit.

The concrete next theoretical search should therefore leave the
entry-regular one-original-site class if retaining the old nonzero
direct V is essential. Alternatively, the canonical regular interaction
already supplies a different, completed Q-activation mechanism. These
are different physical definitions, not interchangeable normalizations
of one measured signal.

## Source and work record

- `7e46c74ce149d5a0a06d1085eb36eebb1bbe6bdb`, execution `branch_only`:
  `notes/local-pair-two-insertion-algebra.md`; canonical regular kernel,
  all single closures, finite-network endpoint theorem and Q activation.
- `2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb`, execution `branch_only`:
  `notes/regular-pair-interaction-result.md` and
  `notes/regular-pair-activation-original-u.md`; completed mixed response,
  full root/slope functional and finite-counterterm dependence.
- `notes/local-pair-size-response-predictions.md`, sections 1--2: actual
  U prefactor and its intrinsic thermal-quotient expression. Its bounded
  occupation-tangent size predictions are not invalidated by this theorem.

This delivery contains definitions and proof only. It performs no test,
enumeration, simulation, root solve, previous-score recalculation or
server operation, and changes no source or coupling in the completed runs.
