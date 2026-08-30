# P250: an exact covariant adaptive-support escape

Date: 2026-08-30. Status: **exact finite-torus witness; pilot interface frozen; no production**.

The fixed-site obstruction is not the end of the physical route.  A minimal state-dependent rule exists that is parameter-free, respects the square-torus geometry, exchanges correctly under the black-NN/white-matching involution, and gives a nonzero order channel already at `L=3`.

## The partial morphisms

A typed marked state is

```text
(omega, h, a_D, a_J, c),
```

where `h` declares the primary NN/matching hand, `a_D` is occupied, `a_J` is vacant, and `c` is a landing mark.

- `D` searches the occupied `h`-graph essential component containing `a_D`.  A candidate must be a non-anchor site which, after deletion, belongs to the vacant complementary component of `a_J`.  Delete the **unique** such candidate closest to `c`.
- `J` searches the vacant complementary-graph essential component containing `a_J`.  A candidate must be a non-anchor site which, after addition, belongs to the occupied component of `a_D`.  Occupy the **unique** such candidate closest to `c`.

The two anchors are protected, so their occupied/vacant types persist through both branches.  This extra transfer condition makes `J` a genuine connector into the black marked component and `D` its exact complement-colour cut operation.

Distance is exact minimum squared Euclidean distance on the axis torus.  If nearest sites tie, the morphism is undefined.  This is important: no coordinate lexicographic convention is used to break a geometric symmetry.

The typed involution is

```text
I(omega,h,a_D,a_J,c)=(C omega,1-h,a_J,a_D,c).
```

It exchanges `D` and `J` site-for-site.  Exhaustive `L=3` checking covers 82,944 defined and undefined support pairs with zero failures.

## Minimal witness

On the `L=3` square torus, occupy the vertical column

```text
omega={(0,0),(0,1),(0,2)}
```

and mark

```text
a_D=(0,0),  a_J=(1,1),  c=(0,1).
```

All four source components are essential.  The supports are

```text
D0          : (0,1)
J0          : (2,1)
J after D   : (0,1)
D after J   : (0,1)
```

Thus deleting the landing site changes the later join support from `(2,1)` to `(0,1)`: the adaptive `J` repairs the cut, whereas `J` first adds the other connector and `D` still removes the landing.  With the projective-leg response at `c`,

```text
L_D=0, L_J=1, L_DJ=1, L_JD=0,
R_plus=0, R_minus=1.
```

So both ordered final fields and the typed response distinguish the two paths.

The same vertical-column witness works at `L=4`.  Applying every translation and all eight square dihedral maps gives 72 exact covariance checks at `L=3` and 128 at `L=4`; every support transforms with the geometry and every witness retains `R_minus=1`.

## Exhaustive `L=3` result

Among states for which all four partial operations are defined, there are exactly 8,136 marked rectangles.  Every one has different ordered final fields and nonzero projective-leg order response:

```text
R_minus=1 : 4,320
R_minus=2 : 3,816
```

The symmetric channel is independently nontrivial:

```text
R_plus=-2 :   360
R_plus=-1 : 3,600
R_plus=0  : 3,168
R_plus=1  :   504
R_plus=2  :   504
```

For each rectangle, its complement+hand-exchanged partner satisfies

```text
R_minus(I state)= R_minus(state),
R_plus(I state) =-R_plus(state).
```

This makes the order channel typed-even and the connected symmetric channel typed-odd for this convention.

## Frozen minimal pilot interface

The exact witness is strong enough to freeze an interface, but not to start stochastic production.  A pilot row must retain:

- geometry, hand, `a_D,a_J,c`, replica ID and base-field digest;
- `D0`, `J0`, `J_after_D`, and `D_after_J` support sites;
- for every selected support: colour, graph, component ID/size, ambient rank, primitive basis, landing distance, and nearest-minimizer count;
- `L_D,L_J,L_DJ,L_JD,R_plus,R_minus`;
- the complement+hand-exchanged paired-row ID.

The nearest-minimizer count is part of the scientific contract: tied rows are excluded, not silently oriented.

This result proves only that one explicit finite-torus physical morphism escapes the fixed-site no-go.  It does not identify a unique continuum operation, prove path memory in existing archives, or imply a Jordan field.

## Reproduction

```sh
python3 scripts/p250_adaptive_support_witness.py
python3 tests/test_p250_adaptive_support_witness.py
```

Output: `results/exact-p250-adaptive-support-witness/latest.json`.
