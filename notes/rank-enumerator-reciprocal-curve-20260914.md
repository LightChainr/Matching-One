# The bivariate topological-rank enumerator and its reciprocal matching curve

Date: 2026-09-14. Exact finite algebraic addendum to draft #773 / #636 / #775.

## 1. One polynomial contains occupation fugacity and topological rank

For a finite N-site torus and primal graph G define the sector fugacity polynomials

    Z_j^G(z) = sum_{S: r_G(S)=j} z^|S|,   j=0,1,2,

and the bivariate rank enumerator

    Q_G(z,u) = Z_0^G(z) + Z_1^G(z) u + Z_2^G(z) u^2.      (1)

Here `z=p/(1-p)` is the ordinary site fugacity and `u` is a topological-rank fugacity. The physical Bernoulli rank-source partition for `X=r-1` differs only by the harmless common factor

    (1+z)^(-N) u^(-1) Q_G(z,u).                           (2)

Thus every finite rank-source quantity used in #773 is encoded by one polynomial which is degree N in z and degree 2 in u.

## 2. Exact reciprocal matching duality

Configuration complement gives

    r_G(S) + r_Ghat(V\S) = 2.

Changing variables `T=V\S` in (1) yields the exact polynomial reciprocity

    Q_G(z,u)
      = z^N u^2 Q_Ghat(z^-1,u^-1).                       (3)

Equivalently, coefficientwise,

    Z_j^G(z) = z^N Z_(2-j)^Ghat(z^-1).                   (4)

For a self-matching model the Newton diagram of Q is centrally symmetric about `(N/2,1)`, and the zero variety is invariant under

    (z,u) -> (z^-1,u^-1).                                (5)

This is the polynomial form of every rank/complement parity relation used elsewhere in the program.

## 3. Canonical rank coordinates are exactly Vieta coordinates of the source roots

For fixed positive real z, the two topological-source zeros are the roots of

    Z_2 u^2 + Z_1 u + Z_0 = 0.                            (6)

Define

    b = 1/2 log(Z_0/Z_2),
    d = log[Z_1/sqrt(Z_0 Z_2)].                           (7)

Writing

    u=e^b v

reduces (6) to the universal centered quadratic

    v^2 + e^d v + 1 = 0.                                 (8)

Hence:

    product(u_+,u_-)=e^(2b),
    v_+ v_-=1.                                           (9)

The two previously introduced canonical coordinates therefore have a direct algebraic meaning:

- b is the logarithmic center of the source-zero pair;
- d controls their separation/type and is invariant under a topological source translation.

The discriminant wall is simply

    d = log 2                                             (10)

or equivalently `Z_1^2=4 Z_0 Z_2`. For `d<log2` the centered roots lie on the unit circle; for `d>log2` they are reciprocal negative reals. The critical modular source-zero collision in the Arguin/Pinson control is precisely the crossing of this algebraic wall.

## 4. Balance is one algebraic section of the same curve

The matching balance condition is

    Z_2(z)=Z_0(z)

and therefore

    b(z)=0.                                               (11)

A topological source `sX` does not alter Q; it merely evaluates the same source curve after translating `log u` by s. Rebalancing under a source is therefore the statement `b(z)=s` derived earlier.

Thus finite matching roots, source tilts and source-zero collisions are not separate polynomial problems. They are different sections of (1).

## 5. Matching parity in canonical coordinates follows immediately

From (4),

    b_G(z) = -b_Ghat(z^-1),
    d_G(z) =  d_Ghat(z^-1).                              (12)

For a self-matching model, with `t=log z`,

    b(t) is exactly odd,
    d(t) is exactly even.                                (13)

This is the fugacity-polynomial explanation of the triangular-site derivative parity controls.

## 6. Long-cylinder tropicalisation

At fixed width and long longitudinal size m, suppose the three normalized rank sectors have exponential rates

    P_0 ~ exp(-m alpha_0),
    P_1 ~ 1,
    P_2 ~ exp(-m alpha_2).

Then the Vieta coordinates obey

    b/m -> (alpha_2-alpha_0)/2,
    d/m -> (alpha_0+alpha_2)/2.                           (14)

Therefore

    alpha_0 = lim (d-b)/m,
    alpha_2 = lim (d+b)/m.                               (15)

The extensive-source three-region phase diagram is simply the tropical/max-plus limit of the three monomials in Q. The width-two open/closed spectral crossing is one explicit local realization of the same geometry.

## 7. Low-fugacity Newton edge for the square torus

For an honest LxL square torus, deterministic geometry gives the leading sector terms

    Z_0(z)=1+...,
    Z_1(z)=2L z^L+...,
    Z_2(z)=L^2 z^(2L-1)+....                             (16)

The first nonzero winding support is a full row or column; a minimum rank-two support is one full row union one full column, sharing one site.

Hence as z->0,

    b(z)= -log L -(L-1/2) log z + o(1),
    d(z)= log 2 + (1/2) log z + o(1),
    c(z)=e^d/2 ~ sqrt(z).                                (17)

Two consequences:

1. the second topological rank costs only `L-1` additional occupied sites after the first, because the two minimum cycles share a site;
2. the centered source-zero angle tends to pi/2 as z->0, showing why the general fixed-source zero-free strip `|Im s|<pi/2` is globally sharp when the whole p range is allowed.

Equation (16) is also a Newton-polygon encoding of the two birth onsets.

## 8. A useful negative result: ordinary real-rootedness fails immediately

It would be convenient if each sector polynomial `Z_j(z)` were real-rooted or belonged to a standard stable-polynomial class. Direct exact-coefficient controls show this is false already at L=3.

For example, square L=3 has

    Z_0(z)=1+9z+36z^2+78z^3+90z^4+45z^5,

which has nonreal conjugate zeros. The same phenomenon occurs broadly in the L=3,4 square/triangular sector polynomials.

Therefore Newton/Lee--Yang analysis of Q cannot assume one-variable sector real-rootedness, interlacing or a strong-Rayleigh shortcut. The exact MLR result for `K|rank2` versus `K|rank0` comes from monotone subset coupling, not polynomial real-rootedness.

## 9. Algebraic research direction

The zero variety

    Q_G(z,u)=0                                             (18)

is a finite "topological spectral curve" with exact reciprocal involution (3). Potentially useful objects are:

- its discriminant in u, `D(z)=Z_1^2-4Z_0Z_2`;
- source-zero branch points `D(z)=0`;
- tropical limits in long cylinders;
- finite-size scaling of the curve near `(p_c,s=0)` after canonical centering.

A continuum near-critical equation of state would amount to a scaling limit of this algebraic curve, not merely convergence of one root.

## 10. Boundaries

- "spectral curve" here is algebraic shorthand; no transfer-operator spectrum or algebraic-geometry novelty claim is implied.
- Sector polynomial non-real-rootedness does not prevent other specialized total-positivity properties; none are claimed.
- The low-fugacity square onset (16) uses honest LxL geometry and should not be transplanted unchanged to degenerate quotients or other lattices.
