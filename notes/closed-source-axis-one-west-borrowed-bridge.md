# The complete one-west class on a rough axis bridge

## Result

A physical one-west contour which borrows its vertical displacement from
the macroscopic directed bridge is parametrically earlier than the closed
L-triomino packet, but still much later than the relaxed west-word bound.
Put

```text
r=m^-1,       c=Lr,       epsilon=Lr^4=L/m^4.               (1)
```

At leading order the unique west edge must be flanked by two vertical
steps of the same sign.  Summing its position and the remaining directed
bridge gives the exact first-order one-particle insertion

```text
delta P_d
 =epsilon [I_(d-2)(2c)+I_(d+2)(2c)].                        (2)
```

Thus the one-west packet has relative weight `epsilon`, not
`beta=L/m^2`.  In the two-boundary noncrossing endpoint ensemble,

```text
J=I0^2-I1^2,
delta J/epsilon
 =2[2 I0 I2-I1(I1+I3)],                                    (3)
```

where every Bessel function has argument `2c`.  Consequently

```text
delta J/J=epsilon C_end(c),
C_end(c)=2[2 I0 I2-I1(I1+I3)]/(I0^2-I1^2),                 (4)

C_end(c)=4-8/c+O(c^-2).                                    (5)
```

The packet is positive and remains `O(epsilon)` conditional on the hard
endpoint.  There is no multiplication by c.  Its common two-boundary bulk
dressing is

```text
C_bulk(c)=4 I2/I0=4-4/c+O(c^-2).                            (6)
```

After removing that common dressing, the genuine Dirichlet distortion is
strictly negative:

```text
C_D(c)=C_end(c)-C_bulk(c)<0,
C_D(c)=-4/c+O(c^-2).                                       (7)
```

Hence the endpoint-specific relative correction is

```text
epsilon C_D(c)=-4/m^3+o(m^-3)              (c->infinity).  (8)
```

This answers the scale question.  The complete leading one-west class is
`beta/m^2=L/m^4` in the endpoint ensemble; its nontrivial Dirichlet part
is `beta/(cm^2)=m^-3`.  It is neither beta nor an alpha-amplified packet.
The sufficient gate contributed by this whole class is `L/m^4->0`.

For the original thermal/source mark, occupation complement pairs the
black and matching-white versions.  At `a=m^-2=0` their centered mark
cancels exactly for every word, before summing (2).  At the physical
two-gas root the first direct source term is therefore

```text
a epsilon M_1(c),                                          (9)
```

where `M_1(c)` is a finite local occupied-corner collar matrix element.
Unlike the shortest L-triomino, the full one-west class contains both
corner placements, so its sign is not fixed by complement alone.  This
single finite collar coefficient is the remaining source datum.  It does
not affect the scale or the negative unmarked Dirichlet result (7).

No `w>=2` contour is included here.

## 1. Exact word count and the missing two powers of m

Fix a based contour word with horizontal homology `(L,0)`, exactly one
west step, n up steps and n down steps.  Its counts are

```text
E=L+1,       W=1,       U=D=n,
T=L+2n+2.                                                   (10)
```

An east neighbour of W would immediately traverse the same cut edge in
the opposite direction.  Both cyclic neighbours of W must therefore be
vertical.  Root the cyclic word at W.  Among the

```text
(T-1)!/[(L+1)!n!^2]
```

rooted words, the exact number with same-sign vertical neighbours is

```text
2n(n-1)(T-3)!/[(L+1)!n!^2].                                (11)
```

Restoring the T possible based origins and dividing by the directed based
count `(L+2n)!/[L!n!^2]` gives the exact fixed-n multiplier

```text
r^2 2n(n-1)(L+2n+2)/[(L+1)(L+2n)].                         (12)
```

The first `r^2` pays for W and its compensating E.  The two forced
vertical neighbours contribute the missing probability `O(r^2)`.  Under
the directed bridge, the Bessel factorial moment is

```text
E[n(n-1)]=c^2 I2(2c)/I0(2c).                               (13)
```

Substitution in (12) yields

```text
Z_(w=1)/Z_directed
 =2 epsilon I2(2c)/I0(2c)+o(epsilon)                       (14)
```

for one boundary, which is (2) at d=0.

Opposite-sign neighbours `U W D` or `D W U` form three sides of a unit
dual square.  With the overwhelmingly likely adjacent E they close that
square and repeat an endpoint/cut edge.  To survive occupied-corner
resolution they need at least one additional adjacent vertical step, and
are `O(epsilon/m)`.  Their formal leading word contribution is a scalar
multiple of `I_d`; such a scalar would cancel from (7) in any event.

The same-sign words have the local forms `U W U` and `D W D`.  Together
with the compensating forward step they are precisely a one-cell corner
flip inserted into a directed macro-edge of vertical displacement `+2`
or `-2`.  This is a planar occupied-corner gadget.  Treating the forward
column as the macro-time makes the decorated network acyclic, so the
ordinary reflection/Lindstrom--Gessel--Viennot determinant applies.  This
is why (2), rather than the relaxed all-word envelope, is the physical
leading insertion.

## 2. Exact endpoint reflection count

Let `P_d=I_d(2c)` be the one-boundary bridge kernel.  At first order in
epsilon, inserting (2) on either one of the two ordered boundaries gives

```text
delta(P0^2-P1^2)
 =2[P0 delta P0-P1 delta P1].                              (15)
```

Since

```text
delta P0/epsilon=2I2,
delta P1/epsilon=I1+I3,
```

equation (15) is exactly (3).  It is positive: it counts the newly
decorated pairs which remain noncrossing.  The same insertion in the
unconstrained two-boundary bulk factor `P0^2` gives (6).

For the sign of their difference, put `q=I1/I0` and
`rho=1-q^2=J/I0^2`.  Adding the harmless scalar `2I_d` turns (2) into
`partial_c^2 I_d`, without changing the normalized rho response.  The
Bessel equation then gives

```text
C_D(c)
 =-2 q^2/(1-q^2)
   [1/c^2-(1/c) partial_c log q].                           (16)
```

The bracket is positive.  Indeed the Turan inequality
`I1^2>I0 I2`, together with `I0-I2=I1/c`, is equivalent to

```text
c(1-q^2)<q,
```

which is precisely `partial_c log q<1/c` through the Riccati equation for
q.  This proves the strict negative sign in (7), not merely its large-c
expansion.  Standard Bessel expansions give (5), (7), and (8).

The mechanism is now explicit.  Most of the one-west weight is a common
positive renormalization of both endpoint and bulk propagation.  The
hard noncrossing wall removes only its diffusion-broadening remainder,
which is smaller by `1/c` and negative.

## 3. Complement and source mark

At zero cloud activity, a decorated rank-one contour and its occupation
complement have the same cut word and the same endpoint determinant.
Their occupied areas are K and `N-K`.  At the centered root `h=1`,

```text
partial_h[h^K+h^(N-K)] after centered normalization =0.    (17)
```

This is wordwise, so neither the reflection subtraction nor the sum over
the insertion position can revive a leading complement-odd mark.

For positive `a=m^-2`, the black and matching-white local gases distinguish
NN and NN+NNN collars.  A same-sign macro-edge changes only a bounded
one-cell corner neighbourhood.  Therefore its failure of (17) is linear
in a with a finite coefficient; it cannot acquire an area or c factor.
Equation (9) is the resulting exact scale definition of `M_1(c)`.

Determining `M_1(c)` requires only the finite list of occupied-corner
collars of the `U W U` and `D W D` gadgets.  Complement symmetry alone
does not determine its sign because the NN and matching collars are not
isomorphic.  This is the minimal unresolved coefficient.  Even if it is
nonzero, its direct source scale is

```text
a epsilon=L/m^6,                                           (18)
```

while the complement-even dressing is `L/m^4`.  Neither can reproduce an
`alpha=L^2/m^3` failure.

## Scientific card

- **Complete leading class:** every physical one-west word whose west
  edge borrows two same-sign vertical bridge steps.
- **Exact scale:** `epsilon=L/m^4=beta/m^2`.
- **Endpoint reflection result:** positive conditional weight
  `epsilon C_end(c)` with `C_end->4`; no endpoint amplification.
- **Dirichlet shape:** after common bulk dressing, the correction is
  strictly negative and equals `-4/m^3+o(m^-3)`.
- **Source mark:** exact complement cancellation at `a=0`; the sole
  remaining datum is the finite collar coefficient `M_1(c)` multiplying
  `a epsilon=L/m^6`.
- **Consequence:** neither the shortest packet nor the complete leading
  borrowed-bridge `w=1` class makes the old alpha gate sharp.  Any alpha
  obstruction must begin beyond this one-west sector.
- **Boundary:** no statement is made about two west steps or interacting
  decorated macro-edges.
