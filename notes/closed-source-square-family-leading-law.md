# Square-family leading law: classify the first two winding layers

**Exact new-size prediction.** Let the first geometry be the axis LxL
torus, L>=5, N=L^2. Let the second be an honest same-area square quotient
with Manhattan period minimum at least L+2. For the existing closed-source
law, lambda=exp(-t), and the original ordered projector denominator
Delta=cos4(theta_axis)-cos4(theta_second),

```text
U(t)/A_N = -(L^2-6L+6)/Delta * lambda^(2L+1)
          + O(lambda^(2L+3)),    A_N=N^(13/8)/2.       (1)
```

In particular the dilated axis/(4,3) pair predicts coefficients
`-46/Delta` at N100, power21, and `-141/Delta` at N225, power31.
The common Delta for these dilations is 1152/625. These follow from a
complete boundary classification and finite-polynomial algebra, without
new histograms, coupling points, fitting or sampling. This companion
continues the [winding-barrier proof](closed-source-winding-barrier.md).

## 1. Minimum layer: only straight stripes

Write gamma=2L-1 and let Z1(h,lambda) be the unnormalized rank-one sum.
The barrier proof gives, for e essential occupied components and c0
zero-image components,

```text
g >= gamma+(e-1)(2L-2)+2c0.
```

At g=gamma there is one occupied component, with exactly two essential
cut-dual boundaries of length L and no other boundary. On the axis torus
a length-L essential boundary must be a straight horizontal or vertical
dual cycle: any other primitive slope costs at least 2L, and a reversal
or transverse step increases the length. The two straight boundaries
enclose a stripe of w full rows, 1<=w<=L-1. Every configuration is counted
once by its orientation, first occupied row and width. Therefore

```text
A(h) := [lambda^gamma] Z1
     = 2L sum_(w=1..L-1) h^(Lw).                    (2)
```

The digital boundary convention remains the occupied-NN one: at an
alternating face the four cut-dual half-edges are paired by turning around
the occupied corners. Every cut edge belongs to one contour; contours
may touch before smoothing but never reuse a cut edge.

## 2. The next layer has exactly two classes

At g=gamma+2, the component bound excludes e>=2 and permits only c0=0
or c0=1. The following two cases exhaust both possibilities.

### Connected case: one partial row attached to a full stripe

For c0=0 the mixed-edge perimeter is 2L+2. Nonaxis winding would cost at
least 4L. The two essential axis contours each have length L modulo2.
There is no budget for a contractible contour, whose minimum length is4,
or for two further essential contours. Their lengths must be L and L+2.
Thus both boundaries cannot be rough at this order: each rough boundary
would separately cost at least two extra edges.

The length-L contour is straight. It forces a full occupied row on one
side and a full vacant row on the other. In particular it cannot pass
straight through an alternating dual vertex, where the stipulated pairing
turns. This already excludes a connected configuration with no full row.

Orient the other contour to have period (L,0). If its east, west, north
and south step counts are E,W,Nv,Sv, then

```text
E-W=L,   Nv=Sv,   E+W+Nv+Sv=L+2,
so W+Nv=1.
```

The option W=1,Nv=0 is impossible: a walk confined to the single cyclic
row with one reversed edge must reuse a cut edge. Hence it has one up
step, one down step, no west step, and L east steps. Their positions
differ by a cyclic distance a in 1..L-1; distance0 or L would again reuse
the vertical edge. The contour is exactly one one-row excursion.

The region between this contour and the straight contour is consequently
w full rows plus a contiguous length-a piece of an adjacent row. Its
parameters satisfy 1<=w<=L-2 and 1<=a<=L-1. The full vacant row supplied
by the straight boundary forbids w=L-1. The representation is unique:
the full occupied rows, partial row, attachment side and cyclic segment
are recovered from the configuration. It cannot be duplicated by the
other orientation, since it has a full vacant row and no full column.

Counting two orientations, L stripe placements, two attachment sides
and L segment starts gives

```text
B_rough(h) = 4L^2 sum_(w=1..L-2) sum_(a=1..L-1) h^(Lw+a).  (3)
```

This contour argument also rules out a hidden alternating-corner family:
the length budget permits only the stated transverse up/down pair.

### Disconnected case: a straight stripe plus one isolated occupied site

For c0=1 all component lower bounds must be saturated. The essential
component therefore has two straight length-L boundaries and is a full
stripe; the other component has mixed-edge boundary4. A finite lifted
NN component with boundary4 is a singleton: two or more cells have at
least six outgoing edges, as follows already from their coordinate
projections. No loop or multi-site appendage has the available budget.

Outside a stripe of width w, exclude its two adjacent rows to prevent
NN attachment. There remain L-w-2 rows, each with L allowed positions.
Thus 1<=w<=L-3 and

```text
B_iso(h) = 2L^2 sum_(w=1..L-3) (L-w-2) h^(Lw+1).    (4)
```

The essential component and the isolated component identify the
representation uniquely. Combining the two cases proves

```text
Z1_axis = lambda^gamma [A(h)+lambda^2 B(h)+O(lambda^4)],
B=B_rough+B_iso.                                    (5)
```

Rank-one g is odd, so there is no intermediate layer. The O term is a
finite polynomial remainder near any fixed positive h.

## 3. Only the isolated-site layer breaks the central occupation symmetry

Both A and B_rough satisfy `P(h)=h^N P(1/h)`. For A this pairs w with L-w;
for B_rough it pairs `(w,a)` with `(L-w-1,L-a)`. Thus
`A'(1)=N A(1)/2` and `B_rough'(1)=N B_rough(1)/2`.
The isolated-site class lacks its complement at this cost: complementing
the isolated occupied island makes a vacant hole in a connected stripe,
reducing C_B by one and increasing g by two.

The finite sums above give the sufficient symbolic values

```text
A0 := A(1) = 2L(L-1),
A'(1) = (N/2) A0,
[d_h^2 (A/(1+h^N))]_(h=1) = -L^4(L^2-1)/6,

B_rough(1) = 4L^2(L-2)(L-1),
B_rough'(1) = (N/2) B_rough(1),
B_iso(1) = L^2(L-2)(L-3),
B_iso'(1) = B_iso(1)*(L^2-L+3)/3,

B'(1)/2-N B(1)/4
  = -L^2(L-2)(L-3)(L^2+2L-6)/12.                   (6)
```

The rough-interface counts are genuine next-layer states, but their
centered derivative cancels. They must not be mistaken for the surviving
asymmetric term in the global leading coefficient.

## 4. Normalization and pooled-root movement supply the other two terms

The only g=2 configurations on either geometry are single occupied sites.
Indeed g>=2C_B+r for a proper occupied set; equality2 forces r=0,C_B=1
and boundary4. Its zero-image lift is finite, and boundary4 forces one
site. A proper rank-two configuration costs at least4, and both rank-one
barriers exceed2. Hence, on each geometry separately,

```text
Z = Z0+Nh lambda^2+O(lambda^4),      Z0=1+h^N,
H_q = h^N-1-Nh lambda^2+O(lambda^4).
```

The zero of the mean of separately normalized q means therefore obeys

```text
h0 = 1+lambda^2+O(lambda^4),     Q_h(h0)=N/2+O(lambda^2).  (7)
```

The positive-coupling simple-root theorem selects this same branch;
there is no fitted root choice. The rank-one probability on the axis is

```text
P1 = lambda^gamma [A/Z0
      +lambda^2 (B/Z0-Nh A/Z0^2)+O(lambda^4)].       (8)
```

Since `[d_h(A/Z0)]_(h=1)=0`, its raw lambda^gamma thermal coefficient
vanishes. With E=1-P1, insert the root movement (7) into (8) to obtain

```text
[lambda^(gamma+2)] E_axis,h(h0)
 = -[d_h^2(A/Z0)]_1 - [d_h(B/Z0)]_1
   + [d_h(Nh A/Z0^2)]_1
 = L^4(L^2-1)/6 - [B'(1)/2-N B(1)/4]
   + N A0(2-N)/8
 = -L^2(L^2-6L+6)/2.                              (9)
```

These are respectively the displaced leading stripe term, the next
rank-one layer, and the common one-particle normalization correction.
None can be dropped while retaining the original root/slope observable.

The companion barrier is at least gamma+4, so its E_h contributes nothing
at the order displayed in (9). Dividing by the ordered Delta and pooled
slope N/2 proves (1), with remainder O(lambda^(gamma+4)). The denominator's
O(lambda^2) correction only enters that remainder.

## 5. New-scale consequence and boundary

For the scaled Gaussian pairs `(5k,0)` and `(4k,3k)`, L=5k and the
companion shortest Manhattan period is 7k>=L+2. Thus

| Area | Pair generators | U/A_N leading term |
|---|---|---|
| 25 | 5, 4+3i | `-lambda^11/Delta` |
| 100 | 10, 8+6i | `-46 lambda^21/Delta` |
| 225 | 15, 12+9i | `-141 lambda^31/Delta` |

For Delta>0 and every integer L>=5 the coefficient is strictly negative:
the original U approaches zero from below at sufficiently strong coupling.
The source barrier itself was 2L-1; the observable begins two powers later
because the minimal stripes are centered at half occupation. The first
nonzero coefficient records both that exact cancellation and the
island-versus-hole asymmetry of the next layer.

This is the fixed-volume, t-to-infinity law of the already specified
microscopic model. It does not interchange the size/coupling limits,
estimate a universal exponent, locate a finite-coupling zero or minimum,
or replace a fresh larger-size numerical readout by independent evidence.
