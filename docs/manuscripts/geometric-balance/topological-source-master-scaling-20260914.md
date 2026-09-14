# A master topological-source scaling function for charge, neutral count, modulus and lattice anisotropy

Date: 2026-09-14

Status: unifying scaling conjecture.  Its finite source coordinates are already exact on the lattice; existence/universality of the continuum scaling function is not proved for square-site Matching One.

## 1. The exact finite generating object already exists

For one honest torus, after aggregating rank-one slope for the moment, define

```text
P0 = Pr(rank=0),
P1 = Pr(rank=1),
P2 = Pr(rank=2),
H_p(z) = E[z^K | rank=1],
```

where `K>=1` is the number of parallel essential components in the neutral/rank-one sector.

Introduce a bounded rank charge `h` and neutral fugacity `z`:

```text
boxed:
Z_L(p;h,z)
 = P0(p) e^-h
   + P1(p) H_p(z)
   + P2(p) e^h.
```

At `z=1`,

```text
partial_h Z|_0 = P2-P0 = M,
```

and the matching root is the zero of this first charge response.

The log-odds coordinate is

```text
vartheta = log(P2/P0).
```

The exact finite common-window decomposition `(chi,vartheta,H)` is simply another coordinate system on the same generating object.

For bond/FK embedded graphs, the new Krushkal observation gives an exact surface-polynomial realization of `h`; the affine-TL seam gives `z=alpha^2/Q` in rank one.  Thus the two source directions already have independent microscopic dictionaries.

## 2. Near-critical continuum target

For a torus of fixed modulus `tau`, let the linear size be `L` and introduce the thermal coordinate

```text
lambda = a_t (p-pc) L^(3/4).
```

The natural continuum target is a topology-resolved scaling function

```text
boxed:
Z_L(p;h,z,tau)
 -> Z_*(lambda;h,z,tau)
```

at fixed `(lambda,h,z,tau)`, after the ordinary normalization appropriate to probabilities.

This one object contains:

```text
rank law                  : z=1, derivatives in h;
neutral essential count   : h=0, derivatives in z within rank one;
matching root             : zero of partial_h Z at h=0;
critical modular geometry : lambda=0, vary tau;
near-critical crossover   : vary lambda;
Potts/TL source dictionary: continue Q and realize h,z by seams/projectors.
```

The main missing continuum engine in #782 is precisely a way to calculate this object (or its sector components) away from `lambda=0`.

## 3. Matching symmetry should act on the source, not only on one observable

At the continuum matched pair, the natural involution is schematically

```text
lambda -> -lambda,
h      -> -h,
z      -> z,
```

with the rank-one projective slope transported by complement.

Thus the strongest form of the continuum symmetry would be

```text
Z_*(lambda;h,z,tau)
 = Z_*(-lambda;-h,z,tau)
```

up to the explicitly declared graph/normalization convention.

Consequences include the familiar critical equality `P0=P2` and the odd/even organization of charge versus neutral variables.

For square SITE this remains a universality/scaling statement, not an exact finite same-graph identity.

## 4. The first lattice correction is naturally a degenerate `x=21/4` shell

For a dimensionless torus source response, the sector-odd correction that produces an `L^-4` root shift has correction exponent

```text
omega_odd = 13/4,
```

because the thermal derivative scales as `L^(3/4)`.

The `x=21/4` resonance note shows that the leading shell can contain at least two angular channels of the same total dimension:

```text
spin 0 : scalar eight-arm-like component,
spin 4 : thermal-family level-four component.
```

Therefore the lattice expansion of the charge source should be written as a **solution space**, not a single ray:

```text
Z_L
 = Z_*
 + L^(-13/4) [
       g0 Z_(21/4,0)(lambda;h,z,tau)
     + g4 Re(e^{i4 theta_lat} Z_(21/4,4)(lambda;h,z,tau))
   ]
 + ... .
```

Here `theta_lat` is the orientation of the microscopic square lattice relative to the chosen cylinder/torus cycle.  Reflection reduces the visible spin-four dependence to a cosine in the current real geometries.

A common sector-even `x≈4` irrelevant field then dresses these coefficients at an additional `L^-2`, producing the observed square `4,6,8,...` root ladder without requiring a new odd operator at every power.

## 5. Root shift from the master expansion

Write the continuum charge log-odds at `z=1` as

```text
vartheta_*(lambda,tau).
```

At the continuum critical point,

```text
vartheta_*(0,tau)=0,
partial_lambda vartheta_*(0,tau) != 0
```

for a nondegenerate thermal coordinate.

Let the `L^-13/4` correction to `vartheta` be

```text
V_odd(tau,theta_lat)
 = b0 F0(tau)
   + b4 Re[e^{i4theta_lat} F4(tau)].
```

Then solving `vartheta_L=0` gives

```text
lambda_root
 = -L^(-13/4)
   V_odd / [partial_lambda vartheta_*(0,tau)]
   + ...,
```

and hence

```text
boxed:
p_root-pc
 = -L^-4
   V_odd /
   [a_t partial_lambda vartheta_*(0,tau)]
   + ... .
```

This formula contains both the scalar and spin-four `x=21/4` amplitudes and makes clear why exponent four alone cannot identify the mechanism.

## 6. Semi-infinite cylinder is a large-aspect boundary condition

For a rectangular torus with aspect `rho=m/w`, the finite charge fugacity obeys at large `rho`

```text
vartheta_(w,m)(p)
 = m Theta_w(p) + O(1),
```

where `Theta_w` is the safe-sector free-energy difference.

At criticality the spin-four cylinder result gives

```text
Theta_w(pc,theta)
 ~ B4 cos(4theta) w^(-17/4)
```

(up to a possible scalar `B0` component now explicitly projected).

Therefore the fixed-aspect torus scaling function for the `L^-13/4` charge correction must have a large-aspect asymptotic whose leading growth is linear in `rho` and whose coefficient matches the semi-infinite-cylinder amplitudes.

This provides a concrete bridge from #585 modular tomography to the much sharper #771 oblique cylinder data:

> candidate torus solution spaces should reproduce the known large-aspect angular projector, not merely fit finite moduli.

## 7. Equal-circumference angular projection extends away from `pc`

The `(4,3),n=2` versus axis `w=10` pair has the same physical circumference.  Once the state spaces are built, evaluating several nearby `p` values is much cheaper than opening more directions.

A stronger test than only `(E,D,p_root)` is to predeclare two or three thermal coordinates and form the angular projector

```text
P4(p)
 = [Theta_axis(p)-Theta_43(p)]/[1-cos4theta_43],

P0(p)
 = [Theta_43(p)-cos4theta_43 Theta_axis(p)]/[1-cos4theta_43].
```

After the appropriate `w^(17/4)` rescaling, the conjecture predicts that these approach two distinct near-critical scaling functions

```text
Psi4(lambda),
Psi0(lambda),
```

rather than merely two constants at `lambda=0`.

This is a high-information extension because it asks whether the angular decomposition is an **operator scaling function**, not just an accidental critical-point fit.

Stop after this pair unless an angular-orthogonal component is actually resolved.

## 8. Fixed-subcritical tail as a topological trans-series

On the subcritical side, define

```text
s = w kappa(p).
```

Near criticality `s` is a function of `lambda` with leading behavior

```text
s ~ C |lambda|^(4/3)
```

in the massive tail.

The rare-topology witness-grading note proposes that the fixed-subcritical expansion is organized by

```text
epsilon = e^-s.
```

Thus the large-negative-`lambda` boundary of the master scaling object should generically admit a topology-graded expansion of the form

```text
observable(lambda)
 = perturbative/massive background
   + sum_(j>=1) e^(-j s(lambda)) P_j(s,lambda),
```

where `j` counts additional essential witnesses/topological instantons and `P_j` carries polynomial/local prefactors.

This is not asserted to be a convergent trans-series.  It is a structural matching rule between:

```text
near-critical finite-lambda scaling
and
fixed-subcritical rare-topology expansions.
```

It explains why pure fragmentation and the coalescing slow-doublet picture are reliable only in the `s->infinity` tail.

## 9. The near-critical genealogy should become split-merge at finite `s`

The #780 no-merger proof gives a relative correction of one-extra-witness order `e^-s`.  Therefore:

```text
s -> infinity : pure cut/fragmentation process;
s = O(1)      : mergers survive and must be resummed.
```

This suggests that the continuum common-label genealogy at finite thermal `lambda` is not a dressed version of the pure fragmentation generator.  It should be a split-merge/topological hard-core process whose far-subcritical expansion reduces to the known cut process.

The same `s` should control the disappearance of the #800 metastable slow doublet.

This supplies a qualitative dynamical boundary condition for any massive Potts/topological process construction.

## 10. Relation to the Krushkal / affine-TL dictionary

For bond/FK configurations on the torus:

```text
h : Krushkal A=e^h, B=e^-h topological source,
z : affine-TL noncontractible-loop seam z=alpha^2/Q.
```

Thus a massive Potts realization of the master object should not start from an arbitrary twist.  It should seek a toroidal modified trace carrying both source coordinates and satisfying the critical homology boundary condition.

The current missing representation-theory step is to realize the Krushkal rank source locally/trace-theoretically in the same periodic algebra that already realizes `alpha`.

## 11. What this reorganizes

This single object puts the current programmes into one hierarchy:

```text
#768 : operator content of the first h-odd lattice correction;
#585 : tau dependence / modular solution space of the same correction;
#767 : lambda dependence and rare-to-critical crossover;
#780 : process interpretation of the large-negative-lambda rare-topology tail;
#782 : massive engine capable of computing Z_*(lambda;h,z,tau);
#800 : one-witness spectral tunnelling in the rare-topology tail.
```

They need not all be active compute tasks.  The point is that a result in one direction now has explicit boundary conditions for the others.

## 12. Minimal falsification programme

1. **Equal-length angular scaling function:** axis `w=10` versus `(4,3),n=2` at a few predeclared thermal points.  If the projected spin-four/scalar decomposition fails immediately away from `pc`, the one-operator-shell picture is too simple.

2. **Large-aspect modular boundary:** any candidate #585 solution basis must reproduce the cylinder `cos4theta` amplitude as `rho->infinity`.

3. **Mass-clock crossover:** compare Palm score speed with independent mass-slope information; failure means `s=w kappa` is not the sufficient clock at the claimed precision.

4. **Finite-s genealogy:** once a controlled sequence with `w kappa=O(1)` exists, measure one merger/split statistic rather than more fixed-subcritical kernel moments.

## 13. Claim boundary

Exact finite ingredients:

- the `(h,z)` source decomposition;
- matching rank/source identities;
- Krushkal finite FK realization of `h`;
- affine-TL finite rank-one realization of `z`.

Conjectural continuum structure:

- existence/universality of `Z_*(lambda;h,z,tau)` for square SITE;
- the two-component `x=21/4` correction shell;
- the topological `e^{-w kappa}` expansion and split-merge finite-`s` process;
- the large-aspect matching between modular and cylinder amplitudes.

The value of this ansatz is that each conjecture now has a different observable direction; failure can localize which bridge breaks.