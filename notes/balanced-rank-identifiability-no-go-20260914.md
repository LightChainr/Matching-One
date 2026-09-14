# Balanced aggregate-rank identifiability no-go

Date: 2026-09-14. Exact finite probability geometry. Addendum to draft #773 and interfaces to #585/#746/#275.

## 1. The balance submanifold is one-dimensional

A three-state rank law is a point

    P=(P0,P1,P2),  P0+P1+P2=1.

At a balance point `P0=P2=a`, `P1=1-2a`. Therefore the entire balanced aggregate-rank family is

    B = { (a,1-2a,a) : 0<a<1/2 }.                       (1)

It is one-dimensional.

In canonical coordinates

    b = 1/2 log(P0/P2),
    d = log[P1/sqrt(P0P2)],

this is exactly the line `b=0`. Under the square-root Fisher embedding

    q_j=2 sqrt(P_j),

`b=0` is the intersection of the radius-2 sphere with the plane `q0=q2` through the origin. Hence B is a Fisher geodesic (a great-circle arc).

## 2. Every balanced tangent has the same rank direction

At `P=(a,r,a)`, `r=1-2a`, an arbitrary centered tangent decomposes uniquely as

    delta P
      = u (-1,0,+1) + v (+1,-2,+1).                      (2)

The first vector is matching-odd; the second matching-even.

If a perturbation is accompanied by whatever compensating source/root adjustment is needed to keep `P0=P2`, then `u=0`. Thus every balanced first-order response of the aggregate rank law is proportional to

    (+1,-2,+1).                                          (3)

There is no second balanced aggregate-rank direction to identify.

This is the tangent version of the exact rank ANOVA basis on #773: the odd coordinate is X=rank-1 and the even coordinate is Y=X^2-E X^2. Rebalancing kills the odd coordinate; only Y remains.

## 3. Stronger statement: the full rebalanced path is source-independent up to parameterisation

Let eta be any physical perturbation for which a compensating scalar source/root can be chosen so that the resulting aggregate law remains balanced. Then necessarily

    P_bal(eta)=(a(eta),1-2a(eta),a(eta)).                 (4)

The unparameterised image of every such experiment is a subset of the **same** geodesic B. Consequently, after forgetting the source parameterisation/amplitude, the full aggregate rank law cannot distinguish:

- a generic-Q deformation after applying the exact topological chemical potential `s_Q=-1/2 log Q`;
- a torus-modulus deformation at critical Q=1;
- any other matching-even source rebalanced to P0=P2;
- a microscopic marked perturbation whose only retained readout is aggregate rank.

If two experiments attain the same value of `a`, their entire aggregate rank distribution is identical. No statistic of rank alone can separate them.

## 4. Direct application to generic Q versus modulus

Arguin/Pinson gives, at criticality,

    P2=Q P0.

After the exact source shift `s_Q=-1/2 log Q`, the law is balanced and is described by one scalar `a_Q(tau)`. At Q=1, varying tau also gives a balanced law described by one scalar `a_1(tau)`.

Therefore the continuum maps

    Q -> P_bal(Q;tau0)

and

    tau -> P(Q=1;tau)

both lie on B. Whenever their a-ranges overlap, there are values `(Q,tau)` producing **exactly the same aggregate rank law** after balancing.

At a high-symmetry modulus such as tau=i, the first derivative with respect to aspect/shear may vanish by modular symmetry while the Q-even tangent is nonzero; that can calibrate a derivative. It does not remove the global identifiability no-go: the target observable space remains one-dimensional.

## 5. Consequence for the graded Q program

The continuum decomposition

    db/d log Q = -1/2,
    dd/d log Q = nontrivial

is still valuable:

1. `db/dlogQ=-1/2` is an exact normalization/homology positive control;
2. `dd/dlogQ` is the one available balanced rank-shape response.

But `dd/dlogQ` cannot by itself identify the microscopic duality-even score with a unique continuum field. Any other balanced even perturbation projects onto the same Y/rank-shape direction.

Thus the exact lattice split of #746 should be interpreted as a **typed source decomposition**, while field identification requires an extra observable beyond aggregate rank.

## 6. Consequence for map-resolved torus tomography

This no-go explains why map/connectivity-resolved observables are not optional if the scientific question is source identification. Aggregate homology rank supplies only one balanced scalar at fixed modulus.

A candidate basis should therefore contain at least one extra coordinate which is not a function of `(P0,P1,P2)` alone, for example:

- a primitive homology direction rather than summed rank one;
- a spatial/Fourier source response;
- a marked connectivity/map sector;
- a conditional structure factor;
- a source-visible boundary/cluster statistic.

Otherwise enlarging the continuum basis only reparameterises the same one-dimensional balanced rank curve.

## 7. Consequence for original-U

The same logic is a finite-volume pre-filter. After removing normalizer and thermal/root nuisance, if two physical candidates differ only through aggregate balanced rank, they are unidentifiable from that readout no matter how much precision is purchased.

A new acquisition is justified only if it adds a represented coordinate outside the rank sigma-algebra / span of the existing nuisance projection.

## 8. Boundaries

- This is an identifiability statement about the **aggregate three-state rank law**, not about the full microscopic measures.
- Source parameterisation/cost can differ even when the unparameterised rank path is the same.
- At fixed high-symmetry modulus, symmetry may make particular nuisance derivatives vanish; that is a calibration, not an increase in target-space dimension.
