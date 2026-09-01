# Width endpoints overturn the all-positive directed capillary candidate

## New result

In the capillary window

```text
L,m -> infinity,       c=L/m in [0,infinity),                  (1)
```

combining the directed boundary transfer with the exact bulk-singleton/
width-zero-mode sum gives the restricted scaling function

```text
Phi_restricted(c)
 = I0(2c)^2/Delta
   *[1-c^2/6+(2c/3) I1(2c)/I0(2c)].                           (2)
```

Here `I0,I1` are modified Bessel functions, equivalently

```text
I0(2c)=sum_(j>=0) c^(2j)/(j!)^2,
c I1(2c)/I0(2c)
 = [sum_(j>=0) j c^(2j)/(j!)^2]
   /[sum_(j>=0) c^(2j)/(j!)^2].                               (3)
```

The earlier version omitted `-c^2/6` and incorrectly claimed strict
positivity. That claim is withdrawn. The corrected bracket is positive
near zero and negative for sufficiently large c, so this restricted model
predicts a capillary sign change in

```text
Ustar/A_N
 =-L^2 m^(-(2L+1)) Phi_restricted(c)[1+o(1)].                 (4)
```

It matches the rigid theorem because `Phi_restricted(0)=1/Delta`.

This is a real signed transmission result for a controlled capillary
submodel, not yet a theorem for the complete rank-one sector. Completeness
requires the explicit interface assumptions in Section 6. The calculation
now keeps three mechanisms together: partition dressing, bridge-area root
motion and width-endpoint singleton leverage. Multiplying the rigid
coefficient by a positive partition function misses the last two.

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

For bulk widths the two transfers factor as in (10). The earlier version
then discarded widths within `O(1)` of zero or L because their width-mode
mass is only `O(1/L)`. That inference is false for the cancellation-sensitive
response: these endpoints carry an `O(L^2)` occupation leverage. Section 2a
keeps the exact bulk-singleton/straight-width correction. A separate narrow-
width interaction of the fluctuating boundaries is not computed here.

This already strengthens the exact two-row obstruction: the smallest
directed closure predicts `I0(2c)^2`, not merely `cosh(c)` or its square.
The partition factor itself is positive; it is not the complete response
coefficient.

## 2a. Bulk singleton gas and the width endpoint

Put

```text
a=c^2/L^2=m^-2,       h=1+a,       B=1+a h=h+a^2.              (10a)
```

For a straight stripe of width w, the available bulk-singleton area in the
restricted geometry is

```text
M_w=max{L(L-w-2),0}.
```

Thus the width sum and its rank-zero/rank-two normalization contain

```text
S_L(h)=sum_(w=1)^(L-1) h^(Lw) B^(M_w),
Z1/(Z0+Z2) proportional S_L(h)/B^N.                            (10b)
```

A direct endpoint expansion at `a=c^2/L^2` gives

```text
d_h log[S_L(h)/B^N] -> c^2/2-c^4/12.                           (10c)
```

To see where the two terms come from, first replace `M_w` by the untruncated
polynomial `L(L-w-2)` for every w. That sum is geometric, and its normalized
logarithmic response tends to `-c^4/12`. Restoring the physical value
`M_(L-1)=0` instead of the fictitious value `-L` changes just one endpoint
term and contributes `+c^2/2`. The endpoint has only `O(1/L)` width weight,
but its occupation leverage is `O(L^2)`. Thus the endpoint-transfer term is
not optional in an `L^2`-normalized, cancellation-sensitive derivative.

Relative to the rigid response, (10c) changes the straight-width coefficient
from `1` to `1-c^2/6`. This is the exact singleton/width-zero-mode correction
inside the present restricted model. It does not yet include the capillary
collision of two interfaces at narrow widths, which can contribute at the
same order.

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

This positive moment supplies the last term in (2). It does not cancel the
independent endpoint-transfer term in (10c).

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

before capillary dressing. Combining the common bulk-width factor `Z_pair`
with the straight-width singleton response (10c) gives the first term

```text
+(L^4/2) m^(-(2L+1)) I0(2c)^2[1-c^2/6].                     (17)
```

There are asymptotically `2L^2 I0(2c)^2` weighted oriented/translated/
width stripes. Equations (13)-(15), followed by the root displacement
`m^-2`, add

```text
+L^4 m^(-(2L+1)) I0(2c)^2 v(c).                              (18)
```

Therefore the corrected restricted prediction is

```text
P1,axis,h
 =(L^4/2)m^(-(2L+1)) I0(2c)^2
   *[1-c^2/6+2v(c)+o(1)].                                    (19)
```

The tilted rank-one sector remains exponentially later. Sector-odds
alignment gives the positive within-geometry denominator

```text
D_h=L^2/2[1+o(1)].                                            (20)
```

Finally `E=1-P1` and `Delta>0` turn (19) into (2)-(4). Root motion changes
the amplitude through (18), while the width endpoint changes it through
(17). Neither changes the winding exponent; slope normalization divides by
a positive quantity and therefore propagates the corrected sign in (2).

## 5. Sign, small-c and large-c predictions

Define the corrected restricted bracket

```text
H(c)=1-c^2/6+(2c/3) I1(2c)/I0(2c).                            (21)
```

It is continuous and `H(0)=1`. On the other hand `I1(2c)/I0(2c)<1`, so

```text
H(c)<1-c^2/6+2c/3<0       when c>2+sqrt(10).                  (22)
```

Consequently the restricted candidate has at least one finite-c zero. This
argument does not establish uniqueness, and a missing narrow-width boundary
interaction can move or remove the zero of the complete rank-one response.

At small c,

```text
I0(2c)=1+c^2+c^4/4+O(c^6),
c I1(2c)/I0(2c)=c^2-c^4/2+O(c^6),
H(c)=1+c^2/2-c^4/3+O(c^6),
Delta Phi_restricted(c)=1+(5/2)c^2+(13/6)c^4+O(c^6).         (23)
```

Thus the capillary cloud initially strengthens the magnitude of the negative
tail, but the endpoint term eventually turns the restricted coefficient. At
large c,

```text
I0(2c)=exp(2c)/sqrt(4pi c)[1+O(c^-1)],
H(c)=-c^2/6+2c/3+O(1),
Delta Phi_restricted(c)
 =-c exp(4c)/(24pi)[1+O(c^-1)].                               (24)
```

Its log growth is `4L/m+O(log(L/m))`, still negligible relative to the
topological tension `2L log m`. Therefore the logarithmic exponent remains
the rigid value even though the amplitude becomes exponentially dressed
on the capillary scale. Because (4) contains `-Phi_restricted`, the restricted
model predicts positive U beyond its candidate zero.

The formula offers two direct prospective failures:

1. the candidate zero is shifted or removed by the omitted narrow-width
   boundary interaction;
2. the normalized capillary amplitude is inconsistent with (2), after the
   known sector/root factors are removed.

Either failure selects missing interface families rather than a failure
of the minimum winding barrier.

## 6. What is complete and what remains a submodel

**Complete in the present restricted calculation:**

- every integer height displacement between successive forward columns at
  bulk widths;
- exact closure/net-zero winding via (7);
- the two independent bulk-width stripe boundaries;
- the exact bulk-singleton/straight-width zero mode, including endpoint
  occupation leverage in (10c);
- the signed bridge-area variance entering root motion;
- the original E sign and positive within-geometry denominator propagation.

**Assumptions needed to identify (2) with the complete rank-one sector:**

1. contours with horizontal reversals/overhangs are `o(1)` relative to the
   directed transfer at fixed c;
2. multi-essential-component and contractible-contour decorations factor
   into the already removed bulk pressure or are lower order;
3. the isolated-site/island-hole asymmetry responsible for the rigid
   coefficient is local apart from the endpoint term already isolated in
   (10c);
4. the narrow-width collision of the two capillary boundaries either
   vanishes at this order or is computed and added explicitly;
5. the tilted entropy remains below its extra `m^-4k` winding cost;
6. restricted sector odds remain aligned so (20) holds.

The first item has a favorable elementary scale: a horizontal reversal
needs a compensating forward step and vertical detour, costing at least
four additional edges; its O(L^2) placements give `O(L^2/m^4)=O(c^4/L^2)`
at fixed c. Similar finite overhang clusters vanish. This does not yet
exclude all mutually interacting overhang families, so it is evidence for
directed completeness rather than a completed theorem.

Items 2-4 are the genuine proof boundary. The old all-positive candidate has
already failed at the straight-width endpoint: an `O(1/L)` sector can survive
an `L^2`-normalized derivative. A capillary boundary collision or a
nonfactorizing island/hole insertion at narrow width can add a signed term
comparable to (17)-(18), and may shift or remove the restricted zero.

## Scientific card

- **New analytic object:** the exact directed height transfer (5)-(9),
  whose scaling partition is `I0(2c)` per interface.
- **Corrected signed candidate:** equation (2) contains the exact
  straight-width factor `1-c^2/6`, is positive near zero and has at least one
  finite-c candidate zero. The earlier strictly positive formula is
  withdrawn.
- **Mechanisms separated:** positive partition dressing `I0^2`, positive
  bridge-area root motion `2v(c)`, and signed width-endpoint transfer
  `-c^2/6`.
- **Next exact target:** compute the narrow-width capillary boundary
  interaction and determine whether it shifts, removes or reinforces the
  restricted zero.
- **Boundary:** the bulk-width transfer and straight-width singleton sum are
  exact in this family; endpoint capillary coupling is incomplete. This is
  not yet the complete square-lattice rank-one sector or a continuum
  interface field.
