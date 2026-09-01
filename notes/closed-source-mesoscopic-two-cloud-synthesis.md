# The closed-source root has a mesoscopic two-cloud expansion

**Main result.**  The strong-source joint limit does not jump directly
from the bounded-Poisson theorem to an uncontrolled thermodynamic regime.
There is an intermediate window in which the rank-zero phase contains a
diverging dilute cloud of occupied sites, the rank-two phase is still the
full configuration up to a vanishing vacancy cloud, and the original
pooled matching root has a parameter-free three-term expansion.

For the equal-area N25 pair, write

```text
N=25 k^2,       a=m^-2,       lambda=N a,       tau=N a^2.
```

If `m -> infinity` and `tau ->0`, then `lambda` may nevertheless diverge.
Uniformly on the root window, the two restricted phases have reference
partition functions

```text
R0(h)=(1+a h)^N,          R2(h)=(h+a^2)^N.                    (1)
```

Their exact equality point is

```text
h_bar=1+a.                                                       (2)
```

After the first connected corrections are included, the pooled Sstar
root obeys

```text
log h_root
 = m^-2 -(1/2)m^-4 -(2/3)m^-6 +o(m^-6).                       (3)
```

The coefficient at `m^-6` is nonzero because the two local pressures
have different first connected counts: one elementary occupied face per
site on the rank-zero side versus two adjacent-vacancy edges per site on
the rank-two side.  Thus `c2-c0=2-1=1`.

This is a new asymptotic prediction for the already named closed source,
not a fit to the N25 exact table.  It extends the existing
`N/m^2 -> zeta < infinity` theorem to a regime with `N/m^2 -> infinity`.
It is still a joint strong-coupling/large-volume theorem, not a fixed-m
thermodynamic statement.

## 1. Complete elementary theorem in the pure-full window

Use the old one-cloud chart

```text
h=(1-a)^-1 exp(s/N),       |s|<=S.                              (4)
```

The growing contour cutoff can be chosen below the systole because
`N/m^4 ->0`.  Keeping the original `exp[A_N(m)]` factor, nested holes and
the occupied-corner resolution gives, for every fixed M,

```text
Z0=(1+h/m^2)^N [1+O_S(N/m^4+N^-M)],                            (5)

Z2/h^N
 =1+N h^-1m^-4
  +O_S(N/m^6+(N/m^4)^2+N^-M),                                 (6)

Z1/(Z0+Z2)=O_S(N^-M).                                         (7)
```

Equation (5) is a relative estimate even when `Z0` grows like
`exp(N/m^2)`.  It is obtained by normalizing nonsingleton black polymers
against the same singleton gas, rather than paying a second volume
factor.  Equation (6) retains the exact N one-vacancy configurations;
all other full-side outer holes start at boundary length six.

Consequently the pooled two-phase law is logistic after centering, both
geometries are individually balanced at the pooled root, and the
within-geometry thermal denominator remains macroscopic.  The existing
rank-one contour estimate then suppresses original U faster than every
inverse power of N along this path.

The conditional structures are asymmetric and explicit:

```text
K | r=0  = Poisson(N/m^2)+o_TV(1),                             (8)

Pr(full configuration | r=2)
          =1-N/m^4+o(N/m^4).                                  (9)
```

The Poisson statement remains meaningful when its mean diverges; it
implies a Gaussian centered/scaled limit and `K/(N/m^2) ->1`.  At the
pooled root the rank-zero and rank-two sectors each have probability
`1/2+o(1)`.  Hence one half of the measure is a mesoscopic particle cloud
even though the other half is asymptotically one full configuration.

The ratio in (9) also proves that `N/m^4 ->0` is the sharp gate for the
claim of conditional concentration on the single full state: the total
one-hole/full weight ratio is exactly `N/(h m^4)`.

## 2. Why the correct center is `1+m^-2`

The rank-zero identity is

```text
g=2K-2 beta_1,
w0=(a h)^K m^(2 beta_1).                                      (10)
```

Every occupied forest therefore has exactly product-site weight.  On the
full background, an isolated vacancy has exact relative activity
`h^-1 a^2`.  Resumming these two singleton clouds gives (1), and

```text
h+a^2=1+a h   iff   h=1+a.                                    (11)
```

The older center `h_c=(1-a)^-1` differs by

```text
h_bar/h_c=1-a^2,
N log(h_bar/h_c)=-N/m^4+O(N/m^8).                             (12)
```

Thus, before the connected polymers are considered,

```text
s_root=-N/m^4+o(N/m^4),
log h_root=m^-2-(1/2)m^-4+o(m^-4).                            (13)
```

The `-1/2` is the first vacancy-cloud correction.  It would be invisible
in an empty-versus-full approximation.

For a finite exact chart, put `r=exp(u/N)` and solve the ideal phase odds
instead of expanding them:

```text
h=(r-a^2)/(1-r a),
(h+a^2)/(1+a h)=r.                                            (14)
```

For Sstar, `u=0` reduces exactly to (2).  A rank-fugacity factor can be
included by replacing `r` with its declared finite-N shifted value; no
new source coefficient is needed.

## 3. The first connected coefficient is also fixed

At `h_bar`, let `theta0` be the occupied probability in the rank-zero
product gas and `theta2` the vacancy probability in the full-background
product gas.  The first Mayer increments are:

```text
rank zero: N (m^2-1) theta0^4 = N/m^6+o(N/m^6),                (15)

rank two: 2N(m^2-1) theta2^2 =2N/m^6+o(N/m^6).                (16)
```

The first line is the one elementary square cycle at each face.  The
second is the one adjacent-vacancy pair at each of the `2N` unoriented NN
edges.  The subtraction by one removes the product-gas contribution that
is already present in (1).  Overlapping cycles, connected vacancy triples
and all other local connected objects start at `N/m^8`.

Therefore the rank-two minus rank-zero log-pressure correction is

```text
d=N/m^6+o(N/m^6).                                             (17)
```

The ideal log-odds slope tends to one, so in the natural chart

```text
u_root=-N/m^6+o(N/m^6).                                       (18)
```

Combining (12) and (18) yields

```text
s_root=-N/m^4-N/m^6+o(N/m^6),                                 (19)
```

and hence (3).  Resolving the displayed constants requires the explicit
long-contour/restricted-sector error to be `o(N/m^6)`; a merely vanishing
bound only proves the centered sigmoid.  The growing-cutoff construction
can target this rate in the `N/m^4 ->0` window.

## 4. A wider two-cloud window and its exact remaining hypothesis

Full-state purity is stronger than root control.  If the vacancy cloud is
retained rather than discarded, both reference gases in (1) remain valid
when `N/m^4` is bounded or divergent.  Their first unresummed interactions
occur at `N/m^6`, not `N/m^4`.

Accordingly,

```text
N/m^6 ->0                                                       (20)
```

is the natural local interaction gate for a two-cloud sigmoid centered at
`h_bar`.  In that wider window the rank-two vacancy count may itself be
Poisson with diverging mean and then Gaussian after scaling; it is not
called a pure full state.

Equation (20) must be paired with a relative long-contour/rank-one bound
below the same phase partition functions.  The old relaxed contour product
contains `exp[O(N/m^2)]` and does not by itself provide that cancellation
throughout the whole wider window.  Thus the local `N/m^6` resummation is
proved, while its unrestricted oblique global extension remains a precise
two-phase contour problem.  In the smaller `N/m^4 ->0` window, the existing
growing cutoff supplies the required relative bound and the complete
claims (3)-(9) follow.

## 5. New shape prediction and stopping boundary

The leading full-side correction is exactly the same N translated
one-hole configurations in both equal-area geometries.  It cancels from
their restricted rank-two shape ratio:

```text
(Z2_f/h^N)/(Z2_s/h^N)
 =1+O(N/m^6+(N/m^4)^2).                                      (21)
```

Therefore a geometry-dependent one-hole amplitude is absent.  The first
possible local shape difference is delayed beyond the common one-hole
cloud; any larger difference must come from connected length-six holes,
large/essential contours or a failure of the declared geometric scope.

This route is complete at the elementary joint-limit level.  It should
not be extended by another fixed-N source mixture or ordinary m64 Monte
Carlo.  The genuinely new continuation would be one of:

1. prove the relative two-phase contour cancellation needed to promote
   the wider `N/m^6` window to the oblique pooled observer;
2. compute the first geometry-sensitive coefficient after (21);
3. connect the exact root coefficients in (3) to the finite-coupling
   polynomial already available, without fitting them.

The first continuation has now exposed a second scaling variable.  The
exact [two-row capillary subfamily](closed-source-capillary-window.md)
has generating function
`[(1+m^-1)^L+(1-m^-1)^L]/2`, so `L/m=sqrt(N/m^2)` controls winding
roughness.  The bulk root theorem remains valid, while the fixed-L U
amplitude requires capillary resummation whenever `L/m` does not vanish.

The theorem does not identify a continuum H4/Jordan field and does not
settle fixed m.  It does identify the correct mesoscopic objects, the
sharp full-purity gate, the first three root coefficients and the first
shape cancellation.

## Dependency card

- **Observer:** original separately normalized pooled q/E root and U.
- **Source:** the fixed closed source Sstar, with no added descriptor.
- **Geometry:** equal-area `(5k,0)` / `(4k,3k)` Gaussian quotients.
- **Input:** existing all-activity contour inequality and exact local
  source weight; no new stochastic or enumeration block.
- **New results:** diverging rank-zero cloud, pure-full coexistence,
  `h_bar=1+m^-2`, noncancelling `c2-c0=1`, expansion (3), and shape
  cancellation (21).
- **Not an independent evidence vote:** all statements are analytic
  consequences of the same fixed model.

Detailed proofs are split into the
[rank-zero gas](closed-source-mesoscopic-black-gas.md),
[rank-two/full phase](closed-source-mesoscopic-full-phase.md), and
[root chart and power gates](closed-source-mesoscopic-root-chart.md).
The separate [capillary window](closed-source-capillary-window.md) explains
why the root limit and the fixed-shell U amplitude have different
uniformity gates.
