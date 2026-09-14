# Critical modular rank-source equation from Pinson/Arguin homology probabilities

Date: 2026-09-14. Continuum positive-control addendum for draft #773 / issue #585. This is a critical-torus calculation, not a near-critical theorem or a lattice field identification.

## 1. Q=1 gives the complete three-state rank law at arbitrary torus modulus

Arguin's torus FK homology formulas (hep-th/0111193, Eqs. (7),(9),(10),(13)) specialize at Q=1 to percolation/Pinson. The homology subgroups are exactly the three aggregate rank classes used in Matching One:

    rank 0 : trivial subgroup {0},
    rank 1 : primitive cyclic subgroup <(a,b)>, summed over primitive directions,
    rank 2 : cross subgroup Z x Z.

At Q=1, the exact continuum duality weight in Arguin Eq. (3) gives

    Z_rank2(tau) = Z_rank0(tau)                            (1)

for every torus modulus tau. Thus the continuum critical rank balance is not special to the square torus: `P0=P2` at every tau.

For Q=1,

    g/4 = 2/3,
    e0  = 2/3.

Dropping a common bosonic prefactor, put

    W_{m,n}(tau)
      = exp[-pi*(2/3)*(m^2 tau_I^2+(n-m tau_R)^2)/tau_I].

Then

    2 Z0 = 2 Z2
         = sum_{m,n} W_{m,n} (-1)^gcd(m,n),                (2)

with gcd(0,0)=0, while the total rank-one weight is

    Z1 = sum_{(m,n) != (0,0)} W_{m,n}
         [ cos(2 pi gcd(m,n)/3) - (-1)^gcd(m,n) ].         (3)

After dividing by `2Z0+Z1`, (2)--(3) give the full critical rank law `(P0,P1,P2)` at arbitrary tau. The implemented Gaussian sums reproduce the square-torus Pinson/Newman--Ziff value already used in #773.

## 2. Every critical modulus supplies an exact topological-source closure

At criticality `b=1/2 log(P0/P2)=0`. Therefore the topological-source partition function is

    Z_tau(s) = P1(tau) + 2 P0(tau) cosh s.                 (4)

The entire critical source response is controlled by the single modular scalar

    c_*(tau) = P1(tau)/(2 P0(tau)).                        (5)

The two source zeros satisfy

    cosh s = -c_*(tau).                                    (6)

Thus:

- `c_*<1`: the nearest pair is purely imaginary, `s=+- i arccos(-c_*)`;
- `c_*=1`: a double static source zero occurs at `s=i pi` modulo `2 pi i`;
- `c_*>1`: the two principal zeros split to `s=+- arcosh(c_*) + i pi`.

This is a static bounded-rank polynomial statement. The double zero is **not** a Jordan block or a thermodynamic phase transition.

## 3. A universal aspect-ratio exceptional point

For rectangular tori `tau=i r`, modular S symmetry gives `c_*(r)=c_*(1/r)`. Solving

    P0(i r)=1/4  <=>  P1(i r)=1/2  <=>  c_*(i r)=1          (7)

with the continuum Gaussian sum gives

    r_* = 1.7878293526796570376290825113182806657619938...
    1/r_* = 0.5593374996898710341376722150641096248339622... (8)

At this universal aspect the critical rank law is exactly

    (P0,P1,P2)=(1/4,1/2,1/4),                              (9)

and the two topological-source zeros collide at `s=i pi`.

For the rectangular controls:

| r | P0=P2 | P1 | c_* | source-zero regime |
|---:|---:|---:|---:|---|
| 1/4 | 0.0710816 | 0.857837 | 6.03417 | real-split + i pi |
| 1/2 | 0.226872 | 0.546257 | 1.20389 | real-split + i pi |
| 1 | 0.309526 | 0.380947 | 0.615372 | purely imaginary |
| 2 | 0.226872 | 0.546257 | 1.20389 | real-split + i pi |
| 4 | 0.0710816 | 0.857837 | 6.03417 | real-split + i pi |

The source-zero topology therefore changes twice under aspect deformation, at `r_*` and `1/r_*`, while the physical percolation system remains at the same critical point. This is an exact warning against interpreting a static zero collision as operator non-semisimplicity.

## 4. Modular positive control for the rank-simplex program

The rank-0 and rank-2 probabilities are modular invariant, while primitive rank-one sectors are permuted by modular transformations before summation. The aggregate `(P0,P1,P2)` is therefore modular invariant. The control script evaluates a generic complex tau and its `T:tau->tau+1` and `S:tau->-1/tau` images and finds agreement at the stored working precision.

This produces a ready-made continuum target for the rank-simplex quantities of #773:

    a_*(tau)=P0=P2,
    d_*(tau)=log[P1/P0],
    rank entropy H_*(tau),
    rank character chi_*(tau)=P0 omega^-1 + P1 + P2 omega,
    topological-source zero locations.

No torus one-point field identification is needed. These are topology-sector partition-function observables already supplied by the Pinson/Arguin continuum theory.

## 5. Interfaces to current tasks

**#585.** Use `(P0,P1,P2)` or equivalently `(a_*,c_*)` across modulus as a mandatory positive control for map-resolved modular tomography. Any candidate implementation unable to reproduce this known homology-sector function should not be trusted on more elaborate matching-odd sources.

**#636.** The source `sX` is a closure-functional deformation only:

    Z_s = e^-s Z_0D + Z_1D + e^s Z_2D.

Arguin's continuum homology weights give its exact critical modular target without changing the local transfer operator.

**#776.** The triangular self-matching finite controls should approach the same `a_*(tau),c_*(tau)` at corresponding torus modulus. This tests topology and modular geometry before any near-critical extension.

## 6. New conjecture: critical modular data as boundary condition for the near-critical equation of state

The near-critical rank-source function proposed in #773,

    Z_tau(lambda,s)
      = Pi0(lambda;tau)e^-s + Pi1(lambda;tau) + Pi2(lambda;tau)e^s,

must satisfy at lambda=0

    Pi_j(0;tau)=P_j^Pinson/Arguin(tau).                     (10)

Together with matching symmetry `Pi2(lambda;tau)=Pi0(-lambda;tau)`, the exact modular curve (2)--(3) is a nontrivial boundary condition for any proposed massive/near-critical continuation.

A useful research target is therefore not an isolated off-critical exponent but the deformation of the known modular scalar `c_*(tau)` into

    c(lambda;tau)=Pi1(lambda;tau)/(2 sqrt(Pi0 Pi2)).       (11)

At lambda=0 this is exactly known for all tau; its lambda derivatives encode the new near-critical information.

## 7. Boundaries

- Arguin/Pinson concerns critical continuum FK/percolation homology probabilities; it does not provide the near-critical lambda dependence.
- The numerical implementation truncates absolutely convergent Gaussian sums; it is high-precision diagnostic arithmetic, not an interval certificate.
- Static source-zero collisions do not imply LCFT Jordan structure.
- No new p_c estimate, Monte Carlo, or original-U identification is made.
