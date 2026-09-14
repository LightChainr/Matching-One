# Rank-softmax geometry: one two-dimensional picture for balance roots and the rank-one plateau

Date: 2026-09-14. Exact finite coordinate geometry for draft PR #773; the asymptotic interpretations are labelled separately.

For any interior three-rank law define

    b = (1/2) log(P0/P2),
    d = log[P1/sqrt(P0 P2)] = log(2c).                     (1)

Then the probabilities are exactly the three-state softmax

    P0 = e^b    / [e^b+e^d+e^(-b)],
    P1 = e^d    / [e^b+e^d+e^(-b)],
    P2 = e^(-b) / [e^b+e^d+e^(-b)].                       (2)

Thus the whole rank law is one point in the `(b,d)` plane.

## 1. All pairwise sector odds are linear

The three log-odds are

    log(P0/P2) =  2b,
    log(P0/P1) =  b-d,
    log(P1/P2) =  b+d.                                    (3)

Hence the three pairwise-equality sets are the straight lines

    rank0 = rank2:   b=0,
    rank0 = rank1:   b=d,
    rank1 = rank2:   b=-d.                                (4)

They meet at `(b,d)=(0,0)`, where all three ranks have probability 1/3.

The ordinary matching balance root uses only the first line `b=0`. It places **no constraint at all** on d. At balance,

    d = log(P1/P0)=log(2c).                                (5)

Therefore two ensembles can have the same perfectly balanced rank-0/rank-2 odds and completely different rank-one dominance.

This is the finite categorical reason a balance-root theorem and a full-law theorem are different objects.

## 2. The physical Bernoulli model is a curve in this plane

For a fixed finite geometry, homogeneous occupation p traces a one-dimensional curve

    p -> (b(p),d(p)).                                      (6)

The existing conditional thermodynamic identities give, in logit z,

    b_z = -(1/2) g1,
    d_z = h1,                                              (7)

where

    g1 = E[K|2]-E[K|0],
    h1 = E[K|1]-(E[K|0]+E[K|2])/2.

Strict monotonicity gives `b_z<0`, so the curve always crosses each vertical b value at most once. The finite matching root is the intersection with the vertical axis.

The two **static sector-equality points**, when they exist, are the intersections with the diagonal rays

    b=d   (P0=P1),
    b=-d  (P1=P2).                                        (8)

These are not being identified with a particular coupling construction's random birth-time medians. They are exact one-time marginal landmarks of the rank law.

## 3. Matching symmetry is mirror symmetry of the physical curve

Primal/matching complement sends

    b -> -b,
    d -> d.                                                (9)

A self-matching model at its symmetric thermal centre therefore has a physical curve invariant under reflection across the d axis. At the centre its tangent is horizontal because

    d_z=h1=0.                                              (10)

This is exact for the triangular-site control at p=1/2. For square site the primal and matching microscopic curves are mirror partners; a universal continuum limit would identify their common limiting shape after the thermal coordinate is aligned.

The odd part of the self-normalized curve is precisely the finite failure of this mirror symmetry within one microscopic lattice description.

## 4. Root versus rank-one plateau

At a balance point b=0:

    P0=P2=1/[2+e^d],
    P1=e^d/[2+e^d].                                       (11)

Thus d has a direct interpretation:

- `d>>1`: the balance root lies deep inside a rank-one-dominated plateau; P0=P2 are both rare;
- `d=0`: all three ranks are equally likely;
- `d<<0`: rank one is suppressed relative to each extreme sector.

This is why equality of the two extreme-sector free energies does not determine the full rank law.

In a fixed-width long cylinder, rank one dominates at the open/closed crossing and d grows extensively with cylinder length. On a fixed-aspect critical two-dimensional torus, d tends a finite topology-dependent constant. These are different limit orders, already visible in the exact width-two bridge and the finite rank-source zero analysis.

### Local approximation to sector separation

Near the balance point, if d varies slowly compared with b, the p-distance between the two sector-equality crossings is approximately

    Delta p_(01,12)
      ~ 2 |d(0)| / |b_p(0)|.                              (12)

This is only a local approximation. The exact distances are obtained by solving `b(p)=+d(p)` and `b(p)=-d(p)`.

Equation (12) makes the two controlling mechanisms explicit:

- d measures the rank-one plateau width in **intrinsic topological-odds units**;
- `|b_p|` converts those units back to the microscopic thermal coordinate.

A balance-root result controls the zero of b. A law-concentration result must additionally control how rapidly the physical curve crosses the other sector boundaries.

## 5. Transfer-spectrum interpretation at fixed width

Suppose for a fixed width and long longitudinal size m the three sector masses have asymptotic forms

    P_j(m) ~ A_j lambda_j^m / Z_total(m).                  (13)

Then common normalization cancels from b,d and

    b/m -> (1/2) log(lambda_0/lambda_2),                   (14)

    d/m -> log lambda_1
           -(1/2)log(lambda_0 lambda_2).                  (15)

So:

- b is the open/closed sector free-energy difference used by cylinder root methods;
- d is the rank-one free-energy excess relative to the geometric mean of the two extreme sectors.

At the cylinder crossing `lambda_0=lambda_2`, b/m=0 but d/m need not vanish. This is the spectral version of the root-versus-law distinction.

Closure prefactors contribute only O(1) terms to b,d at fixed width; they matter for finite-length root shifts but not the leading per-length limits.

## 6. Topological sources are translations in this phase plane

Couple sources `sX+tX^2`. Up to a common normalization the three sector weights become

    (e^(b-s+t), e^d, e^(-b+s+t)).                          (16)

After removing the common factor e^t from the extreme sectors, the canonical coordinates transform as

    b -> b-s,
    d -> d-t.                                              (17)

Thus `(b,d)` are affine source coordinates for the full interior rank simplex. The three equality lines (4) are the exact categorical phase boundaries under these finite closure sources.

The ordinary topological source (`t=0`) translates horizontally; the even extreme-vs-middle source translates vertically.

This geometric statement concerns a bounded closure variable. It is not a bulk phase transition or a local transfer Jordan structure.

## 7. Continuum target

For a near-critical scaling limit at torus modulus tau, the proposed universal curve can be written directly as

    lambda -> (b_*(lambda;tau),d_*(lambda;tau)),           (18)

with matching symmetry

    b_*(-lambda;tau)=-b_*(lambda;tau),
    d_*(-lambda;tau)= d_*(lambda;tau).                    (19)

The self-normalized equation `d=D_*(b;tau)` eliminates lambda and its thermal metric entirely. Its odd finite-size residual is the parameterization-free matching-asymmetry diagnostic described in the companion note.

This two-dimensional representation is a useful organizational object for #613/#739 (root versus law), #636 (closure weights), #768 (matching-odd corrections), #775 (exact rank-sector polynomials) and #776 (self-matching control). It does not replace any of their proof obligations.
