# Exact canonical coordinates for the three-rank simplex

Date: 2026-09-14. Additive exact note for draft PR #773. No scaling hypothesis is needed in Sections 1--3.

Let

    P_j = P(r=j), j=0,1,2,
    X=r-1 in {-1,0,1},

with all three probabilities positive. Define

    b = (1/2) log(P_0/P_2),
    c = P_1 / [2 sqrt(P_0 P_2)] > 0.                       (1)

These two real numbers give a global coordinate chart on the interior of the rank-probability simplex.

## 1. Exact inverse map

Write `t=sqrt(P_0P_2)`. Then

    P_0=t e^b,
    P_2=t e^(-b),
    P_1=2ct.

Normalization gives `2t(c+cosh b)=1`, hence

    P_0 = e^b / [2(c+cosh b)],
    P_1 = c   / [c+cosh b],
    P_2 = e^(-b) / [2(c+cosh b)].                          (2)

Thus `b` measures rank-0/rank-2 imbalance while `c` measures the rank-1 weight relative to the geometric mean of the two extreme sectors.

The ordinary matching balance condition is simply

    b=0.                                                    (3)

## 2. The complete rank-source partition function

For a topological source `s` coupled to `X`,

    Z_X(s)=P_0e^(-s)+P_1+P_2e^s.

Using (2),

    Z_X(s) = [c+cosh(s-b)]/[c+cosh b].                     (4)

Hence every complex source zero is determined by

    cosh(s-b)=-c.                                          (5)

The two coordinates have a direct zero-geometry meaning:

- `b` is the horizontal centre of the two zero families;
- `c` controls their separation from that centre.

If `0<c<1`,

    s=b +/- i arccos(-c) + 2 pi i n.

If `c=1`, the two families coalesce into double zeros at `s=b+(2n+1)i pi`. If `c>1`,

    s=b +/- arcosh(c) + (2n+1)i pi.

This repackages the general rank-source zero formula without losing information: `(b,c)` reconstructs the full rank law.

## 3. Topological-source action is a pure translation

Tilt the law by `exp(sX)` and renormalize. Then

    P_0 -> P_0 e^(-s)/Z_X(s),
    P_1 -> P_1/Z_X(s),
    P_2 -> P_2 e^s/Z_X(s).

Substitution into (1) gives exactly

    b -> b-s,
    c -> c.                                                (6)

Thus a pure topological chemical potential translates only the matching-odd coordinate `b`; the even shape coordinate `c` is an exact invariant of that source. The source needed to rebalance a configuration ensemble is therefore `s=b`.

This is a useful implementation control: a code path advertised as coupling only to `X=r-1` must leave `c` unchanged at every finite size and every real source.

## 4. Matching/complement parity

Under primal/matching complement exchange,

    (P_0,P_1,P_2) -> (P_2,P_1,P_0),

so

    b -> -b,
    c -> c.                                                (7)

Therefore the proposed near-critical continuum equation of state can be reduced from the triplet `Pi_j(lambda;tau)` to two scalar functions

    b(lambda;tau), c(lambda;tau),

with the symmetry prediction

    b(-lambda;tau)=-b(lambda;tau),
    c(-lambda;tau)= c(lambda;tau).                         (8)

This makes the roles of the two pieces explicit:

- `b` contains the matching-odd thermal/topological response and determines the balance root;
- `c` is an independent matching-even rank-shape observable and determines the imaginary zero angle at balance.

At critical balance, `b=0` and

    c_* = (1-2a_*)/(2a_*) = 1/(2a_*)-1,                   (9)

where `a_*=P_0=P_2`. For the square-torus Pinson/Newman--Ziff value `a_*=0.3095262754298313...`,

    c_* = 0.61537174608411717427...,
    theta_* = arccos(-c_*) = 2.233653823652404... .       (10)

## 5. Thermal derivatives in canonical coordinates

With logit thermal coordinate `z`, sector thermodynamics gives

    b_z = (1/2)(E[K|0]-E[K|2]) = -(1/2) Delta E[K],       (11)

and

    (log c)_z
      = E[K|1] - (1/2)(E[K|0]+E[K|2]).                    (12)

Equation (11) is exactly the inverse susceptibility appearing in the topological-source root slope:

    dz_*/ds = 1/b_z = -2/Delta E[K].                       (13)

Under a smooth matching-symmetric near-critical limit, `b(lambda)` is odd and `c(lambda)` even. Thus `(log c)_z` should have no leading `L^(3/4)` contribution at the symmetric continuum centre; a surviving large term is a direct diagnostic of finite matching-odd corrections or an incorrect thermal-centre convention.

The second derivative of `b` is the rank-sector conditional variance difference,

    b_zz = -(1/2) Delta Var(K),                             (14)

and the higher derivatives alternate similarly. This coordinate system therefore packages the conditional-cumulant hierarchy into one odd function `b` plus one even function `c`.

## 6. Scaling use and claim boundary

A natural future universal target is not a single fitted exponent but the pair of continuum functions `(b(lambda;tau),c(lambda;tau))` after an independently calibrated thermal coordinate. Dimensionless Taylor combinations such as

    b'''(0)/b'(0)^3,
    c''(0)/b'(0)^2

are invariant under an overall rescaling of `lambda` and can be compared across lattices or map-resolved continuum constructions.

This note does not prove the square-site near-critical scaling limit, identify a CFT field, or modify the original-U source. It supplies exact finite coordinates and typed targets for those later questions.
