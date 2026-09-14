# Fisher area versus canonical persistence: an exact geometric bridge

Date: 2026-09-14. Exact inequality and finite controls. Addendum to draft #773.

## 1. The two global quantities come from the same self-normalized rank curve

Use the canonical coordinate

    b=1/2 log(P0/P2)

and the random canonical persistence interval `[B2,B1]` from the two rank births. Its mean length is

    E G = E(B1-B2) = integral_R P1(b) db.                 (1)

The Fisher-sphere area between the physical rank-law curve and the `P1=0` geodesic is

    A_F = 2 integral_R sqrt(P1(b))/cosh(b) db.            (2)

Both are invariant under every reparameterisation of the microscopic occupation coordinate.

## 2. Exact persistence lower bound from Fisher area

Cauchy--Schwarz and

    integral_R sech^2(b) db = 2

give

    (A_F/2)^2
      <= [integral P1(b) db] [integral sech^2(b) db]
      = 2 E G.

Therefore

    E G >= A_F^2/8.                                      (3)

This is an exact inequality for every finite rank curve with positive sector probabilities; no scaling or percolation-specific input is used.

Equality would require

    P1(b) proportional to sech^2(b)                      (4)

almost everywhere. The actual percolation curves need not satisfy (4); the saturation ratio below measures their distance from this extremal shape in the corresponding weighted L2 sense.

## 3. A new dimensionless shape invariant

Define

    eta_pers = A_F^2/[8 E G],  0<eta_pers<=1.             (5)

Unlike a root shift, amplitude or raw birth gap, this scalar depends only on the unparameterised self-normalized rank equation of state.

High-precision controls from the exact L=3,4 rank-sector polynomials are:

| lattice | L | E G | A_F^2/8 | eta_pers |
|---|---:|---:|---:|---:|
| square |3|1.111923978620|1.023345704706|0.920337833|
| square |4|1.173080389706|1.073809574831|0.915375949|
| triangular |3|1.149474505170|1.058181140146|0.920578173|
| triangular |4|1.180362367040|1.084524300787|0.918806234|

The cross-lattice spread is already small at these sizes. This is a positive universality signal, not an asymptotic estimate.

**Conjecture FP1.** At fixed torus modulus, `eta_pers` converges to a lattice-independent continuum value together with the full self-normalized rank curve.

This is a global shape test: two curves can agree at `b=0` and in several Taylor coefficients but disagree in eta_pers.

## 4. Spherical isoperimetric constraint

The physical rank curve plus the `P1=0` geodesic edge forms a closed loop on the radius-2 Fisher sphere. If its enclosed area is the chosen `A_F` branch, the spherical isoperimetric inequality gives

    (L_F + pi)^2 >= A_F [4 pi - A_F/4],                  (6)

where `L_F` is the Fisher length of the physical curve and `pi` is the geodesic distance between the pure rank-0 and pure rank-2 endpoints on the radius-2 sphere.

Equation (6) is not expected to be close to equality; it is a consistency region for the pair `(L_F,A_F)` independent of lattice details.

## 5. Process interpretation

Equation (3) connects a one-time information-geometric observable to an actual two-birth persistence statistic. It does not reconstruct the birth copula: `E G` is itself determined by the one-time interval-coverage function `P1(b)`.

Higher moments of G still require multi-time persistence information. Thus the hierarchy is:

    rank curve -> E G and Fisher geometry,
    two-time kernel -> E G^2 and copula shape,
    full paired births -> entire G law.                   (7)

This makes eta_pers a useful bridge between #775/#776 (one-time exact rank curves) and #778 (paired-birth process).

## 6. Boundaries

- Numerical controls use exact coefficient polynomials plus high-precision quadrature, not interval-certified integrals.
- The small-size cross-lattice agreement is not a proof of universality.
- eta_pers is a property of the aggregate rank curve, not a continuum field identifier.
