# The shortest occupied-corner hairpin packet

## Result

The relaxed west-step bound does not describe the shortest physical
square-lattice packet.  Fix one straight horizontal cut edge from
`(0,0)` to `(1,0)`.  A local simple replacement which contains exactly
one west step, stays on one side of the untouched straight contour, and
has no repeated cut edge has length seven.  Up to translation along the
edge, the four words are

```text
upper: U U E E D W D,       U W U E E D D,
lower: D D E E U W U,       D W D E E U U.                 (1)
```

Each word, closed by the deleted straight edge, is the boundary of an
L-triomino.  The two upper words are its left/right placements and the
lower words are their vertical reflections.  There is no length-three or
length-five physical packet.  Consequently a shortest hairpin changes

```text
Delta Bmix=6,       Delta C_B=Delta r=0,       Delta g=6.    (2)
```

Its bare local activity is `m^-6`, and its translated bulk activity is
`O(L/m^6)`, not the relaxed `L/m^2` activity.  At a width-one endpoint the
packet directed into the thin strip is forbidden, while its outward mate
is legal.  Therefore the unmarked endpoint contact coefficient is
strictly negative.  It does not cancel.

The **original thermal/source mark does cancel at leading order**, for a
different reason.  The two complementary thin endpoints retain a black
L-triomino protrusion and a matching-white L-triomino indentation.  With
the actual NN-black/NN+NNN-white collar convention their exact local
increments are

```text
                 Delta K     Delta H     Delta M
black protrusion     +3           0          -6
white indentation    -3          -8          +1,            (3)
```

where H counts sites eligible for an isolated matching-white hole and M
counts sites eligible for an isolated NN-black singleton.  Put

```text
a=m^-2,       C=h+a^2,       A=1+ah.
```

After the common carrier factor is removed, one lateral placement has

```text
w_plus(h)=h^3 A^-6,
w_minus(h)=h^5 C^-8 A.                                    (4)
```

At the two-gas root `h=1+a`, `A=C=1+a+a^2`.  Hence the order-`a^0`
thermal derivatives of (4) are `+3` and `-3`: the complementary endpoint
pair cancels exactly.  The first residual is

```text
[d_h(w_plus+w_minus)]_(h=1+a)
 =a(1+a)^2(3a^3+4a^2-4a-8)/(1+a+a^2)^8
 =-8a+O(a^2).                                               (5)
```

Thus the shortest packet has a negative direct source residual, but only
at `O(L/m^8)`.  Its unmarked dressing remains `O(L/m^6)`.  It cannot
saturate the earlier possible `alpha=L^2/m^3` failure scale.  This is an
exact exclusion of the shortest local packet, not yet a theorem for a
west step whose vertical sides are borrowed from a macroscopically rough
bridge.

## 1. Local-word classification

Let a replacement path start at `(0,0)`, end at `(1,0)`, and contain one
west step.  Since the net horizontal displacement is one, it contains two
east steps.  A length-three candidate `E,E,W` necessarily immediately
retraces a cut edge.  Adding one up/down pair gives length five.  The west
edge must have vertical neighbours, but the remaining path then either
closes a unit dual square before reaching the endpoint or meets the
undeformed straight contour at an interior vertex.  Both violate the
resolved-contour rule.

With two up and two down steps, the replacement has length seven.  Add
the deleted edge in the reverse direction.  The result is a simple dual
polygon of perimeter eight lying entirely on one side of that edge.
Because it uses two east, one west, two up and two down steps, the polygon
cannot be a straight triomino or a `2 x 2` square with another edge on the
old contour.  It is an L-triomino.  Its marked edge has two lateral
placements and two choices of side, giving exactly (1).

This argument also explains the occupied-corner role.  The other formal
seven-step self-avoiding walks cross or touch the undeformed straight
collar away from the marked edge.  At the resulting four-cut dual vertex,
turning around occupied corners either separates a contractible contour
or reuses a cut edge; it does not produce a one-carrier replacement.

Equation (2) now follows without a source approximation.  The detour has
six more cut edges than the deleted edge.  Flipping the three enclosed
sites attaches to the existing carrier (or to its matching-white
complement), so neither the number of essential occupied components nor
the ambient rank changes.

## 2. Endpoint sign before applying the source mark

For a bulk-width annulus all four words in (1) are legal at either
boundary.  At gap one, a word lying in the inter-boundary strip intersects
the other resolved contour and is rejected by the noncrossing gate.  The
outward word remains legal.  For one boundary edge the resulting local
two-particle contact insertion is therefore

```text
-2m^-6 w_in |d=1><d=1|,                                   (6)
```

where the factor two is the left/right L-triomino placement and `w_in` is
the appropriate positive weight from (4).  Equation (6) is negative and
nonzero.  Equivalently, the shortest packet lowers endpoint survival
relative to the bulk packet-dressed transfer.  There is no cancellation
between the two lateral words: both remove positive paths.

This is the precise boundary of the scalar-renormalization intuition.
Away from contact, the seven-edge detour begins and ends at the same two
vertices as the deleted edge and merely dresses that local propagation.
At `d=1`, half of this scalar packet is projected out.  That projection is
the negative Dirichlet/contact coefficient.

## 3. Exact source and collar weights

The two dilute local factors have different microscopic neighbourhoods.
An eligible matching-white hole must have all eight NN+NNN neighbours
occupied.  An eligible NN-black singleton must have its four NN
neighbours vacant.  Take the straight black half-plane to be `y<=0`.
The two protrusions are

```text
{(0,1),(0,2),(1,2)},       {(0,1),(0,2),(-1,2)},
```

and the two complementary indentations remove

```text
{(0,0),(0,-1),(1,-1)},     {(0,0),(0,-1),(-1,-1)}.
```

Checking only these sites and their eight-neighbour collars gives (3):

- a black protrusion destroys six exterior singleton positions and
  creates no eight-neighbour interior hole position;
- a white indentation destroys eight interior hole positions and creates
  one exterior singleton position.

The full local gas factor is

```text
h^(K-H) C^H A^M.                                           (7)
```

Substitution of (3) into (7) gives (4).  This is why the indentation has
`h^5`, rather than the incorrect bare factor `h^-3`: removing three
occupied sites simultaneously removes eight factors which had already
contributed a power of h through the hole gas.

At `a=0`, (4) becomes `h^3` and `h^-3`.  Their derivatives at `h=1` are
opposite.  Equation (5) is the exact finite-a remainder at the physical
two-gas root.  Both lateral placements multiply it by two and leave its
sign unchanged for sufficiently small positive a.

There are therefore two distinct statements:

1. the unmarked endpoint contact is negative and survives;
2. its leading original-source projection cancels between the two
   complementary thin endpoints, leaving a negative `O(a)` mark.

The first dresses the already negative endpoint law by relative
`O(L/m^6)`.  The second is a genuinely new direct source channel of order
`O(L/m^8)`.  Neither has the `O(alpha)` scale obtained by dividing the
relaxed west-word weight by a generic endpoint probability.

## Scientific card

- **Exact packet:** the four words (1), equivalently the two lateral
  L-triominoes on either side of a marked straight edge.
- **Exact cut cost:** `Delta g=6`; translated activity `L/m^6`.
- **Endpoint geometry:** the inward half is forbidden, giving a strictly
  negative unmarked Dirichlet/contact insertion.
- **Source mechanism:** black protrusion and matching-white indentation
  have unequal collars but cancel at order `a^0`; the residual is
  `-8a+O(a^2)` per complementary lateral pair.
- **Consequence:** the shortest packet cannot make `alpha=L^2/m^3` sharp.
  A possible alpha crossover, if any, must use a nonlocal west step whose
  vertical sides are supplied by the rough bridge rather than this local
  L-triomino packet.
- **Boundary:** no claim is made here about all one-west contours or their
  resummation.
