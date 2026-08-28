# Exact staggered tangent family and the improved-action target

Status: exact `N=10` tangent calculation plus a frozen `N=130` discovery protocol.
This note constructs an explicit local matching action. It does **not** identify
the resulting tangent with a unique continuum primary.

## Local family

Let `e` denote the degree-8 (`x+y` even) sublattice and `o` the degree-4
sublattice of the C4 checkerboard triangulation. Occupy sites independently with

\[
p_e=\frac12+t+\lambda,\qquad p_o=\frac12+t-\lambda .
\]

The exact probability domain is the diamond

\[
|t+\lambda|\leq\frac12,\qquad |t-\lambda|\leq\frac12.
\]

Because the graph is its own site-matching graph, occupation complement gives
an exact local pair map

\[
F(t,\lambda)=(-t,-\lambda),\qquad J=dF|_0=-I_2,
\]

and `F^2=1` without an OPE assumption. The `t` direction is the ordinary
thermal tangent. The staggered `lambda` direction is a new scalar microscopic
tangent which distinguishes the inequivalent 8/4 coordination sublattices but
preserves C4.

For any wrapping channel define

\[
R^\pm(t,\lambda)=\frac12\{R(t,\lambda)\pm R(-t,-\lambda)\}.
\]

Then `R+` and `R-` are exact pair-exchange eigenobservables. A Taylor monomial
`t^m lambda^n` belongs to `R+` when `m+n` is even and to `R-` when `m+n` is odd.
This is a lattice tangent grading, not a claimed local-CFT/OPE grading.

## Exact minimum-quotient response

`scripts/c4_self_matching_tangent.py` exhausts all `2^10` configurations of
the `(a,b)=(3,1)` quotient and groups masks by the two sublattice occupation
counts. It therefore constructs `R(t,lambda)` as an exact bivariate rational
polynomial, not a finite-difference estimate.

For all five recorded channels (`cross`, `both`, `either`, `direction_0`, and
`direction_1`) the center response, with rows `(R+,R-)` and columns
`(t,lambda)`, is

\[
B=\left.\frac{\partial(R^+,R^-)}{\partial(t,\lambda)}\right|_0
=
\begin{pmatrix}
0&0\\[1mm]
15/8&5/4
\end{pmatrix}.
\]

It obeys the explicit intertwining relation

\[
\operatorname{diag}(1,-1)B=BJ.
\]

Thus thermal and staggered perturbations are both exchange odd at first order;
the response of an exchange-even observable vanishes, while the exchange-odd
response is nonzero. At second order the selection reverses.

Along the critical-center slice `t=0`, the exact odd polynomial is channel
independent:

\[
R^-(0,\lambda)=\frac54\lambda-4\lambda^5
=\lambda\left(\frac54-4\lambda^4\right).
\]

Its real roots are `lambda=0` and

\[
|\lambda|=(5/16)^{1/4}=0.747674\ldots .
\]

The latter lie outside the legal interval `|lambda|<=1/2`. Therefore the
minimum exact quotient gives a sharp negative result:

> The only legal zero of the leading finite-quotient exchange-odd response is
> the self-matching center `lambda*=0`; there is no nonzero odd improved point
> in this one-parameter probability family.

This does not prove that every larger quotient has no extra root. It does make
`lambda*=0` the only exact candidate and turns any claimed nonzero odd root
into a directly falsifiable finite-size phenomenon.

The exchange-even polynomials remain nontrivial. For example,

\[
R^+_{cross}=\frac{11}{32}+\frac58\lambda^2+rac52\lambda^4-10\lambda^6,
\]

\[
R^+_{either}=\frac{21}{32}-\frac58\lambda^2-rac52\lambda^4+10\lambda^6,
\]

while `R+_direction_0=R+_direction_1=1/2` on this quotient. This is the reason
the nontrivial improved-action search should target the even H4 amplitude, not
try to manufacture a second odd root.

## Frozen nontrivial H4 search

For two inequivalent representations of the same `N`, define

\[
A^+_{T4}(N,\lambda)=
N\,\frac{R^+_{(a_1,b_1)}(\lambda)-R^+_{(a_2,b_2)}(\lambda)}
{\cos(4\theta_1)-\cos(4\theta_2)}.
\]

The prefactor `N` tests the leading identity-family `x=4` correction. Exact
exchange makes this an even function of `lambda`. The `N=10` quotient cannot
define this contrast because it has no inequivalent same-norm representation;
the first compatible design is

```text
N=130: (11,3) versus (9,7).
```

Freeze the nonnegative discovery grid

```text
lambda = 0, 1/8, 1/4, 3/8.
```

Use common uniform fields across orientations. For each field at `lambda`,
also score its occupation complement, whose law is `-lambda`; their half-sum
is the antithetic estimator of `R+`. Fit only the first three points to

\[
A^+_{T4}=a_0+a_2z+a_4z^2,\qquad z=\lambda^2,
\]

and reserve `lambda=3/8` as a no-refit lack-of-fit point. A nonzero candidate is
admissible only when the fitted discriminant is nonnegative, a root satisfies
`0<z*<1/4`, the held-out residual passes, and a fresh antithetic run locally
brackets zero around `lambda*=sqrt(z*)`. The decisive replication transports
that frozen root, without refitting, to

```text
N=170: (13,1) versus (11,7).
```

Failure to find an admissible root is a useful result: this staggered tangent
then calibrates matching parity but is not an improved action for the leading
even spin-4 lattice coupling.

## Reproduction

```bash
python3 scripts/c4_self_matching_tangent.py \
  --json results/local-20260828/P44-self-matching-tangent/exact-tangent.json
python3 -m unittest discover -s tests -p 'test_c4_self_matching_tangent.py' -v
```
