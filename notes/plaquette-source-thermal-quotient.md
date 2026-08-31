# The fixed plaquette source modulo the common thermal clock

**Exact mechanism statement.** The single-site part of the full-plaquette
source is a common change of the Bernoulli thermal coordinate. Its response
to the root-comoving, slope-normalized original U vanishes exactly. Therefore
a nonzero V_F4 must come from its centered two-, three- or four-site
interaction part; no sign or continuum-field identity follows from this
algebra.

Work on the ordinary square torus with N sites and N unit faces, with
injective3×3 stencils so all face corners and the pair types below are
distinct. The Gaussian bound N>8 suffices. The source and positive law are

```text
F4(eta)=sum_faces product_(v in face) eta_v,
P_(p,lambda)(eta) proportional to P_p(eta) exp(lambda F4(eta)).
```

F4 is fixed when p is differentiated. All identities are at lambda=0 and
0<p<1. No sampling, old residual readout or fitted response is used here.

## 1. One forced face gives the exact source response

For a fixed translation-invariant configuration observable O, put
`O_bar(p)=E_p[O]`, `O_square(p)=E_p[O|one specified face full]`, and
`delta_O=O_square−O_bar`. Independence and translation invariance give

```text
J_O(p)=Cov_p(O,F4)=N p^4 delta_O(p),
J'_O(p)=N [4p^3 delta_O(p)+p^4 delta'_O(p)].
```

Indeed each face indicator has probability p^4 and the same conditional
expectation. After forcing its four distinct sites, only N−4 independent
sites remain. Thus O_square has the exact Bernstein representation

```text
O_square(p)=sum_(k=0)^(N−4) BinomPMF(k;N−4,p) O_square,k,
O'_square=Cov_p(O,K_outside|face full)/[p(1−p)].
```

The conditional polynomial has degree at most N−4; delta_O need not, since
it subtracts the unconditional degree-at-most-N polynomial. This is a
Bernoulli conditional law, not an unconditional N-site convolution with
four threshold indices shifted afterward. The identities apply separately
to each geometry and then to the fixed pooled/projected q/E combinations.

## 2. The exact centered interaction expansion

Let xi_v=eta_v−p and K=sum_v eta_v. Expanding each face gives

```text
F4 = Np^4 + 4p^3(K−Np) + p^2 Q2 + p Q3 + Q4,
Q2 = 2 sum_NN_edges xi_u xi_v + sum_face_diagonals xi_u xi_v,
Q3 = sum_faces sum_(three corners of face) product xi_v,
Q4 = sum_faces product_(four corners) xi_v.
```

Each site occurs in four faces. Each NN edge occurs in two faces; each of
the two diagonals of a face occurs in that face only. Thus Q2 contains2N
NN edges with multiplicity2 and2N diagonals with multiplicity1:6N pair
incidences, exactly six per face. Q3 has4N right-angle triples, each with
multiplicity1; Q4 hasN plaquettes. These counts use the nonalias scope above,
not a silently simplified small quotient.

Independence makes each Q_j centered and orthogonal to every single-site
function, hence to K. They are the exact Bernoulli interaction orders2,3,4.
For distinct sites S, `Cov(O,product_(v in S) xi_v)` also equals the connected
joint cumulant of O with those site variables; their mutual independence
kills every proper partition not containing all of them with O.

Since `Cov(O,K)=p(1−p) O'_bar`, the single-site source response is

```text
J_O^single = 4p^4(1−p) O'_bar = f(p) O'_bar,
J_O^multi = Cov(O,p^2 Q2+p Q3+Q4),
J_O = J_O^single + J_O^multi.
```

The same f(p) holds for both geometries: there is no extra factor N in this
clock coefficient. The bulk-N factor is already present in the forced-face
formula, consistently with the full sum F4 rather than F4/N.

## 3. The original global-U differential kills every common clock

Write Q(p)=mean(q), Y(p)=P4(E), with the original fixed geometry weights,
and evaluate on the simple pooled root Q(p0)=0. Here Q denotes the matching
mean, not cluster fugacity. Put

```text
A=A_N (independent of p and lambda),  D=Q',  r=Y'/D,
U=A Y'/D,   j_Q=partial_lambda Q,   j_Y=partial_lambda Y.
```

The root moves by `a=−j_Q/D`. Differentiating numerator and denominator on
that same root gives the exact functional

```text
L[j_Q,j_Y] = A/D * [j'_Y−r j'_Q−r' j_Q]
          = A/D * partial_p(j_Y−r j_Q),
r'=(Y''−r Q'')/D.
```

In particular, for **any** differentiable common function f(p),

```text
j_Q=f Q',   j_Y=f Y'  =>  j_Y−r j_Q=0  =>  L[j_Q,j_Y]=0.
```

This includes f' and both root/slope terms; dropping the derivative of r
would destroy the identity. Applying it to the expansion proves
`V_F4=L[J^multi_Q,J^multi_Y]` exactly. If the forced-face route is used,
the same derivative can be evaluated directly as

```text
V_F4 = A N/D * {4p^3(delta_Y−r delta_Q)
                +p^4(delta'_Y−r delta'_Q−r' delta_Q)} at p0.
```

The centered xi_v and their coefficients depend on p. Either differentiate
all of that dependence in J^multi, or freeze the entire decomposition at
p0. The latter single-site term has coefficient4p0^3 and response
`4p0^3 p(1−p) O'`, which is still a common clock and is likewise annihilated.
Mixing the two conventions by freezing xi but differentiating only its
displayed p coefficients is not the derivative of the fixed F4 experiment.

## 4. Why this source is specified, and what a nonzero result would mean

The [checkerboard endpoint map](square-checkerboard-endpoint-homology.md)
proves `Ctot_parent=Ctot_child+F4`: a filled child plaquette represents an
isolated forced parent-sublattice site. Thus the additional action lambda F4
is dictated by decimation of the absolute cluster source. It is a static
local interaction in a positive occupation measure, distinct from the
stopped rank-conditioned lag-one P154 source.

A nonzero V_F4 would establish transmission of this prescribed multisite
interaction to the original global U beyond a common density clock. A zero
would not remove the interaction: L can cancel other source directions too.
Multisite interactions can still couple to occupation polynomials and
thermal descendants; the result would not identify thermal-Q4, a continuum
primary, a field count or the unique source of the unperturbed anomaly.
Nothing here predicts the sign, changes the frozen source or proposes an
additional observable or fit.
