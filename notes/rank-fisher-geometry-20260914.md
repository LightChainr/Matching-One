# Rank-simplex Fisher geometry: intrinsic curvature, information retention, and global near-critical fingerprints

Date: 2026-09-14. Addendum to draft PR #773. Exact finite algebra plus high-precision controls at square/triangular L=3,4. No exponent fit or field identification.

## 1. The three rank probabilities form a constant-curvature statistical manifold

For `P=(P0,P1,P2)` in the positive simplex, use the Fisher--Rao metric

    ds^2 = sum_j dP_j^2/P_j.

The square-root map

    q_j = 2 sqrt(P_j)

is an isometric embedding into the positive octant of the sphere `sum q_j^2=4`. Hence the rank-probability manifold has Gaussian curvature exactly `1/4`. This is model-independent finite probability geometry.

Use the canonical natural coordinates already introduced in #773,

    b = 1/2 log(P0/P2),
    d = log(P1/sqrt(P0 P2)) = log(2c).

Then

    P = (e^b, e^d, e^-b)/(2 cosh b + e^d),

and the Fisher metric is the Hessian of

    A(b,d)=log(2 cosh b + e^d).

At a balance point `b=0`, writing `P0=P2=a` and `P1=r=1-2a`,

    g_bb = 2a,
    g_dd = 2ar,
    g_bd = 0.                                             (1)

Thus the matching-odd coordinate b and matching-even rank-shape coordinate d are exactly Fisher-orthogonal at every finite balance root.

## 2. The physical rank-law curve has an intrinsic geodesic curvature

Eliminate the microscopic occupation coordinate and write the physical curve as

    d = D_L(b).

This is the same parameterisation-free object as `c=C_L(b)`, since `d=log(2c)`.
At `b=0`, put

    s = D'_L(0),
    u = D''_L(0).

The only Christoffel terms needed at the balance line are

    Gamma^b_bd = -r/2,
    Gamma^d_bb = -1/2,
    Gamma^d_dd = (1-2r)/2.

A direct covariant-acceleration calculation gives the signed Fisher geodesic curvature

    kappa_F
      = sqrt[r/(2a)] * [u - 1/2 + s^2/2] / (1+r s^2)^(3/2).   (2)

This is invariant under every reparameterisation of p or the thermal scaling field. In the self-matching case `s=0` exactly; in a matching pair the conjectured continuum symmetry gives `s->0`.

The combination `u-1/2`, rather than u alone, is important: a Fisher geodesic tangent to the b direction itself has `d''=1/2` in these natural coordinates. Therefore the previously observed `d'' about -0.11` is a genuinely large intrinsic bending of the physical topological equation of state, not a coordinate artefact.

Using the exact rank-sector polynomials already committed on #773:

| lattice | L | D'(0) | D''(0) | signed kappa_F |
|---|---:|---:|---:|---:|
| square | 3 | 0.0624022 | -0.125653 | -0.448504 |
| square | 4 | 0.0337935 | -0.112253 | -0.453036 |
| triangular | 3 | 0 exactly | -0.122896 | -0.458651 |
| triangular | 4 | 0 exactly | -0.117056 | -0.461399 |

Two tiny sizes do not establish a limit, but the cross-lattice agreement is substantially tighter than one would expect from raw p-coordinate derivatives. This is a high-value L=5--8 continuation for #775/#776.

## 3. Exact rank ANOVA: there are only two nontrivial rank-information channels

At a balance root define

    X = r_rank - 1 in {-1,0,1},
    Y = X^2 - 2a.

Then

    E X = E Y = 0,
    Cov(X,Y)=0,
    Var X = 2a,
    Var Y = 2ar.                                           (3)

`X,Y` span all centered functions of the three-state rank. Hence for any square-integrable observable A,

    E[A | rank] - E A
      = Cov(A,X)/(2a) X + Cov(A,Y)/(2ar) Y,

and exactly

    Var(E[A|rank])
      = Cov(A,X)^2/(2a) + Cov(A,Y)^2/(2ar).                (4)

There is no hidden third rank channel. The first term is the matching-odd/topological channel; the second is the matching-even/rank-shape channel. This is the finite-probability counterpart of the typed split used in #746, but it requires no continuum identification.

## 4. Thermal information retained by rank

Take A=K, total occupied-site number, and use logit z as the thermal coordinate. Let

    g1 = E[K|r=2]-E[K|r=0],
    h1 = E[K|r=1] - 1/2(E[K|r=0]+E[K|r=2]).

Then the Fisher information about z carried by rank alone is

    I_rank = Var(E[K|rank])
           = a g1^2/2 + 2ar h1^2.                         (5)

The full Bernoulli configuration has

    I_full = Var K = N p(1-p).

Therefore the exact information-retention fraction is

    eta_rank = I_rank/[N p(1-p)].                          (6)

The odd part is exactly the earlier thermal--topological squared Fisher correlation; the even part is the additional rank-shape information. In the self-matching triangular control `h1=0` exactly, so all rank thermal information lies in the odd channel.

If `g1~L^(3/4)` and the matching-even slope is smaller, then

    eta_rank ~ const * L^(-1/2).                           (7)

Thus rank remains an O(1) near-critical topological readout while retaining a vanishing fraction of the microscopic Bernoulli thermal information. This gives a theoretical reason not to expect K-only conditioning to remain efficient asymptotically in the cut-network program.

Small-size controls:

| lattice | L | odd fraction | even fraction | sqrt(L) eta_rank |
|---|---:|---:|---:|---:|
| square |3|0.653924|0.000870|1.13414|
| square |4|0.580388|0.000235|1.16125|
| triangular |3|0.712444|0|1.23399|
| triangular |4|0.615471|0|1.23094|

The tiny even square contribution is itself another matching-asymmetry observable.

## 5. Mutual information conjecture

The complete rank-sector occupation polynomials give the exact mutual information

    I(rank;K).

At the finite balance roots:

| lattice | L=3 | L=4 |
|---|---:|---:|
| square | 0.579831 | 0.446867 |
| triangular | 0.665358 | 0.490754 |

The products `sqrt(L) I(rank;K)` are approximately 1.004, 0.894 (square) and 1.152, 0.982 (triangular).

**Conjecture F1.** Under a rank-conditioned density CLT with common leading variance, `I(rank;K)=Theta(L^-1/2)`. In the weak-separation Gaussian limit its leading term should be `eta_rank/2`; this stronger coefficient statement is not yet supported at L<=4 and is explicitly falsifiable with #775.

## 6. Two global Fisher invariants of the entire near-critical rank law

Because the Fisher metric is intrinsic, the total length of the physical path from p=0 (pure rank 0) to p=1 (pure rank 2),

    L_F = integral sqrt(sum_j (dP_j)^2/P_j),                (8)

is independent of the thermal parameterisation. The endpoints are orthogonal points of the radius-2 Fisher sphere, whose geodesic distance is exactly pi. Hence

    L_F >= pi,                                             (9)

with equality only for the boundary geodesic that never occupies rank 1.

A second invariant is the radius-2 spherical area between the physical rank curve and the `P1=0` great-circle edge. With spherical coordinates

    sqrt(P0)=sin(theta) cos(phi),
    sqrt(P1)=cos(theta),
    sqrt(P2)=sin(theta) sin(phi),

one has `phi=atan(exp(-b))`, and the area is

    A_F = 4 integral sqrt(P1) d phi
        = 2 integral_R sqrt(P1(b))/cosh(b) db.              (10)

It therefore depends only on the self-normalized rank equation of state, not on p.

Executed high-precision quadrature gives:

| lattice | L | Fisher length L_F | L_F-pi | Fisher-sphere area A_F |
|---|---:|---:|---:|---:|
| square |3|4.04065843|0.89906578|2.86125246|
| square |4|4.07381003|0.93221738|2.93095148|
| triangular |3|4.06266217|0.92106952|2.90954449|
| triangular |4|4.07995877|0.93836612|2.94553805|

**Conjecture F2.** At fixed torus modulus these converge to lattice-independent continuum values `(L_F^*,A_F^*)`. Unlike a fitted exponent or p-coordinate amplitude, both are global functionals of the whole near-critical rank law and contain no thermal metric factor.

This is a stronger universality target for #775/#776 than matching one or two Taylor coefficients.

## 7. Boundaries

- The Fisher geometry is exact finite probability theory; the claimed cross-lattice continuum equality is a conjecture.
- Numerical quadrature above uses exact coefficient polynomials but is not interval-certified.
- Neither Fisher curvature nor global length identifies a CFT field.
- `eta_rank -> 0` concerns information about the microscopic thermal parameter in the rank coarse-graining; it does not say the rank variable becomes trivial in the near-critical topology problem.
- No new Monte Carlo or STATUS claim is introduced.
