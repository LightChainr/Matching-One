# A directed capillary transfer gives a signed candidate `Phi(c)`

## New result

In the capillary window

```text
L,m -> infinity,       c=L/m in [0,infinity),                  (1)
```

the smallest closed transfer model which contains every directed
axis-boundary height excursion has an explicit original-U scaling
function

```text
Phi_DSOS(c)
 = I0(2c)^2/Delta * [1+(2c/3) I1(2c)/I0(2c)].                 (2)
```

Here `I0,I1` are modified Bessel functions, equivalently

```text
I0(2c)=sum_(j>=0) c^(2j)/(j!)^2,
c I1(2c)/I0(2c)
 = [sum_(j>=0) j c^(2j)/(j!)^2]
   /[sum_(j>=0) c^(2j)/(j!)^2].                               (3)
```

Every term in (2) is nonnegative and the bracket is strictly positive.
Therefore this model predicts

```text
Ustar/A_N
 =-L^2 m^(-(2L+1)) Phi_DSOS(c)[1+o(1)] <0,                   (4)
```

with no finite-c zero. It matches the rigid theorem because
`Phi_DSOS(0)=1/Delta`.

This is a real signed transmission result for a controlled capillary
submodel, not yet a theorem for the complete rank-one sector. Completeness
requires the explicit interface assumptions in Section 6. The value of
the exercise is that it identifies both the partition dressing and the
root-motion area term; multiplying the rigid coefficient by a positive
partition function alone would miss the second factor in (2).

No configuration, simulation, finite-c fit or source change is used.

## 1. Exact finite transfer for one directed boundary

Represent an axis essential boundary by its integer height after crossing
each of the L horizontal columns. Require one forward horizontal crossing
per column but allow an arbitrary vertical displacement between successive
columns. If `r=m^-1`, the displacement kernel is

```text
T_r(a,b)=r^|a-b|,       a,b in Z.                              (5)
```

Its Fourier eigenvalue is exact:

```text
t_r(theta)=sum_(d in Z) r^|d| exp(i d theta)
          =(1-r^2)/(1-2r cos(theta)+r^2).                     (6)
```

After fixing the irrelevant global height, the closed-loop partition is
the zero-displacement matrix element

```text
Z1,L(r)=[z^0](sum_(d in Z)r^|d|z^d)^L
        =(1/2pi) integral_0^(2pi)t_r(theta)^L dtheta.           (7)
```

For `r=c/L`,

```text
L log t_(c/L)(theta)=2c cos(theta)+O(c^2/L),                   (8)
```

uniformly at fixed c. Hence

```text
Z1,L(c/L) -> (1/2pi) integral exp[2c cos(theta)]dtheta
            =I0(2c).                                         (9)
```

The earlier exact two-row family `cosh(c)` is contained in this transfer:
it keeps only two heights, while (5) permits both signs of every unit
height step. Displacements of absolute size at least two have total rate
`O(L/m^2)=O(c^2/L)` and vanish in the limit, but retaining them in (5)
makes the finite transfer exact before the limit.

## 2. Two stripe boundaries

For a stripe whose width is a positive fraction of L, the two essential
boundaries fluctuate in disjoint O(1)-height neighborhoods with probability
tending to one. Their directed transfers therefore factor:

```text
Z_pair(c)=I0(2c)^2.                                          (10)
```

Widths within O(1) of zero or L occupy only `O(1/L)` of the width zero
mode and do not affect the L^2-normalized response. Boundary collision is
therefore absent at leading order in this declared model. The orientation,
translation and width multiplicities remain those of the rigid stripe;
(10) is their capillary multiplier.

This already strengthens the exact two-row obstruction: the smallest
directed closure predicts `I0(2c)^2`, not merely `cosh(c)` or its square.
All are positive, so partition dressing alone cannot reverse the sign.

## 3. The bridge-area moment is also exact

The root shift differentiates the occupation area, so the partition factor
is insufficient. In the limit behind (9), condition on j upward and j
downward unit jumps. Their weight is `c^(2j)/(j!)^2`. Conditional on j,
the upward and downward jump locations are two independent sets of j
uniform points on the unit circle. If H(x) is the height bridge, its signed
area is

```text
mathcal A=integral_0^1 H(x)dx
         =sum_(down times)t - sum_(up times)t,
E[mathcal A|j]=0,       Var(mathcal A|j)=j/6.                  (11)
```

Thus one interface has

```text
Var(mathcal A)
 =c I1(2c)/[6 I0(2c)].                                       (12)
```

The occupation-area change between two independent stripe boundaries is
L times the difference of their bridge areas. Consequently

```text
Var(delta K)/L^2
 =v(c):=c I1(2c)/[3 I0(2c)].                                 (13)
```

This positive moment is the extra root-motion term in (2).

## 4. Propagation through the actual root and denominator

Let `P1(h)` denote the separately normalized axis rank-one probability.
For a reciprocal rank-one family,

```text
P1,h(1)=0.
```

In log activity `u=log h`, a configuration pair with occupation K and
N-K contributes

```text
d_u^2 [sum h^K/(1+h^N)]_(u=0)
 =1/2 sum [(K-N/2)^2-N^2/4].                                 (14)
```

Adding capillary area fluctuations changes (14) by

```text
(1/2) number_of_weighted_stripes * Var(delta K).               (15)
```

The two-cloud root has `h_root-1=m^-2+o(m^-2)` at the order of the
winding response. The rigid fixed-L cancellation, after root motion,
next-shell asymmetry and ordinary normalization, gives

```text
P1,axis,h
 =+(L^4/2)m^(-(2L+1))[1+o(1)]                                (16)
```

before capillary dressing. Multiplying all local rigid insertions by
`Z_pair` gives the first term

```text
+(L^4/2) m^(-(2L+1)) I0(2c)^2.                               (17)
```

There are asymptotically `2L^2 I0(2c)^2` weighted oriented/translated/
width stripes. Equations (13)-(15), followed by the root displacement
`m^-2`, add

```text
+L^4 m^(-(2L+1)) I0(2c)^2 v(c).                              (18)
```

Therefore the directed-SOS prediction is

```text
P1,axis,h
 =(L^4/2)m^(-(2L+1)) I0(2c)^2[1+2v(c)+o(1)].                 (19)
```

The tilted rank-one sector remains exponentially later. Sector-odds
alignment gives the positive within-geometry denominator

```text
D_h=L^2/2[1+o(1)].                                            (20)
```

Finally `E=1-P1` and `Delta>0` turn (19) into (2)-(4). Root motion changes
the amplitude through (18), but not the winding exponent; slope
normalization divides by a positive quantity and cannot change the sign.

## 5. Sign, small-c and large-c predictions

The series in (3) make positivity elementary:

```text
Phi_DSOS(c)>0       for every c>=0.                            (21)
```

At small c,

```text
I0(2c)=1+c^2+O(c^4),
c I1(2c)/I0(2c)=c^2+O(c^4),
Delta Phi_DSOS(c)=1+(8/3)c^2+O(c^4).                          (22)
```

Thus the capillary cloud strengthens the magnitude of the negative tail
in this model; it does not initially bend it toward a zero. At large c,

```text
I0(2c)=exp(2c)/sqrt(4pi c)[1+O(c^-1)],
Delta Phi_DSOS(c)=exp(4c)/(6pi)[1+O(c^-1)].                    (23)
```

Its log growth is `4L/m+O(log(L/m))`, still negligible relative to the
topological tension `2L log m`. Therefore the logarithmic exponent remains
the rigid value even though the amplitude becomes exponentially dressed
on the capillary scale.

The formula offers two direct prospective failures:

1. a sign change of the fully resummed coefficient at finite c;
2. a normalized capillary amplitude inconsistent with (2), after the
   known sector/root factors are removed.

Either failure selects missing interface families rather than a failure
of the minimum winding barrier.

## 6. What is complete and what remains a submodel

**Complete inside the directed-SOS model:**

- every integer height displacement between successive forward columns;
- exact closure/net-zero winding via (7);
- both independent stripe boundaries and the width zero mode;
- the signed bridge-area variance entering root motion;
- the original E sign and positive within-geometry denominator.

**Assumptions needed to identify (2) with the complete rank-one sector:**

1. contours with horizontal reversals/overhangs are `o(1)` relative to the
   directed transfer at fixed c;
2. multi-essential-component and contractible-contour decorations factor
   into the already removed bulk pressure or are lower order;
3. the isolated-site/island-hole asymmetry responsible for the rigid
   coefficient is local and receives the common factor (10);
4. noncrossing interaction of the two boundaries and narrow widths is
   lower order after the L^2 normalization;
5. the tilted entropy remains below its extra `m^-4k` winding cost;
6. restricted sector odds remain aligned so (20) holds.

The first item has a favorable elementary scale: a horizontal reversal
needs a compensating forward step and vertical detour, costing at least
four additional edges; its O(L^2) placements give `O(L^2/m^4)=O(c^4/L^2)`
at fixed c. Similar finite overhang clusters vanish. This does not yet
exclude all mutually interacting overhang families, so it is evidence for
directed completeness rather than a completed theorem.

Items 2-4 are the genuine proof boundary. In particular, a capillary-
dependent island/hole insertion which does not factor through (10) can add
a signed term comparable to (18). Such a term is the only currently named
route to a zero of the full `Phi(c)`; the positive partition and area terms
in (2) cannot create one.

## Scientific card

- **New analytic object:** the exact directed height transfer (5)-(9),
  whose scaling partition is `I0(2c)` per interface.
- **New signed candidate:** equation (2), strictly positive for all c, so
  the original U remains negative with no capillary zero in this model.
- **Mechanism separated:** ordinary partition dressing `I0^2` and the
  positive root-motion area response `2v(c)`.
- **Prospective discriminator:** a full rank-one resummation sign change
  must come from a nonfactorizing island/hole or coupled-boundary term, not
  from directed capillary roughness itself.
- **Boundary:** this is complete directed-SOS rank-one transfer, not yet the
  complete square-lattice rank-one sector or a continuum interface field.
