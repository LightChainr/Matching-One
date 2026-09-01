# Cloud and source closure under the beta gate

## Result

Let

```text
N=L^2,       a=m^-2,       c=L/m,       beta=L/m^2,
m->infinity,                    beta->0.                     (1)
```

After the common bulk polymer pressures and the already retained straight
two-gas collar are factored out, the complete contractible cloud changes
the noncrossing endpoint transfer and its original thermal mark by

```text
O(ac)+O(Na^2)_common+O(ac/m^2)+O(exp[-kappa ell log m]).     (2)
```

The `Na^2=beta^2` term is a same-area, geometry-independent pressure.  It
is a common normalization/thermal-clock change and disappears from the
angular original-U quotient.  The first genuinely geometry-dependent
remainder is therefore

```text
ac=c/m^2=beta/m=o(1).                                      (3)
```

It comes from the black-singleton eligibility collar of the rough
interface.  Matching-white holes and all connected boundary polymers are
smaller by another factor `m^-2`.  The same bounds survive one thermal
derivative, so (3) closes the cloud/source part of the beta gate.

There is no additional collar state of size `sqrt(c)`.  The height range
of a capillary bridge is `O(sqrt(c))`, but eligibility is local: only a
one-site-thick tubular neighbourhood of the cut contour changes.  Its
cardinality is controlled by the number of vertical cut edges, whose
endpoint expectation is `O(c)` with exponential tails.  Multiplying the
range by an entire row length would count sites which remain locally
eligible and is not the square-lattice collar.

The conclusion is an error bound for the actual closed source.  It does
not introduce an independent dilute-gas observable or assume that the
common pressure itself is small; its leading singleton part may be
extensive and is removed exactly.

## 1. Exact endpoint control of the rough collar

In the continuous directed limit, the two-boundary bulk and endpoint
partitions are

```text
G(c)=I0(2c)^2,
J(c)=I0(2c)^2-I1(2c)^2.                                   (4)
```

Mark every unit vertical edge by replacing c with `c exp(t)`.  The total
vertical variation S then satisfies

```text
E_bulk S=c partial_c log G(c),
E_end  S=c partial_c log J(c).                              (5)
```

The Bessel formulas give, uniformly for `c>=1`,

```text
E_bulk S=4c-1+O(c^-1),
E_end  S=4c-2+O(c^-1),                                    (6)
```

and for c in a compact interval both are bounded.  The same marked
partition at a fixed small positive t gives

```text
E_end exp[tS] <= exp[C_t(1+c)].                             (7)
```

Thus all fixed moments of S are `O((1+c)^k)`.  In particular endpoint
conditioning does not turn `E S=O(c)` into `O(c^2)` by a crude division by
`J/G~1/c`; the determinant can and should be differentiated directly.

For the exact arbitrary-run kernel the marked symbol is obtained by
`r^|d| -> exp(t|d|)r^|d|`.  Its logarithmic derivatives differ from (5)
by `O(1/m+L/m^3)`, so the same `O(1+c)` collar bound holds throughout the
beta window in which the directed kernel is used.

Now compare a rough contour with its straight carrier of the same width.
Every site whose singleton/hole eligibility changes lies within bounded
graph distance of a vertical edge or a turn.  A cut edge has bounded
degree and can be charged only a bounded number of such sites.  Therefore

```text
|Delta collar| <= C(1+S).                                  (8)
```

This is the promised exclusion of a separate `sqrt(c)` state.  The range
may visit `sqrt(c)` rows, but (8) counts the contour tube itself, not full
rows.  Equations (5)-(8) give an endpoint collar expectation `O(1+c)` and
a rough-minus-straight expectation `O(c)`.

## 2. The two local gases and their exact activities

Use the same notation as the complete axis transfer:

```text
z_B=ah,             A=1+z_B,
z_W=a^2/h,          C=h+a^2=h(1+z_W).                       (9)
```

Here `z_B` is the activity of an NN-black singleton in the vacant
background and `z_W` that of an isolated matching-white hole in an
occupied background.  At the pooled capillary root h stays in a fixed
compact subset of `(0,infinity)`, so

```text
z_B=O(a),       z_W=O(a^2).                                (10)
```

The two straight collar rows and their endpoint truncation are already
included exactly in the factors A and C and in the hard kernel J.  Only
the eligibility difference (8) remains.  Changing one eligibility bit
changes the logarithm of its local factor by `O(z_B)` or `O(z_W)`.  Hence

```text
|Delta log Z_collar|_end
 <=C E_end[(1+S)(z_B+z_W)]
 =O(a+ac)+O(a^2+a^2c).                                    (11)
```

The constant-a term is a finite endpoint convention and vanishes with m.
The first unbounded-c term is

```text
ac=beta/m.                                                 (12)
```

No generic `1/J` amplification is needed: (11) is already an expectation
inside the normalized endpoint ensemble.  Equivalently, differentiating
`log J` shows that the leading `O(ac)` vertical-length dressing is common
to bulk and endpoint; their difference is only `O(a)`.  The weaker bound
(12) is sufficient and remains uniform without using that cancellation.

## 3. Connected contractible polymers

The common cloud is not restricted to independent singleton sites.  It
has a standard local positive polymer expansion at large m.

### Vacant-background black polymers

A connected NN-black component of size at least two has mixed perimeter
at least six.  Including its component reward, its activity begins at
`O(m^-4)`.  The number of rooted lattice animals of size k is exponential
in k with a fixed lattice constant, while every additional exposed edge
supplies another inverse power of m.  Consequently the pressure beyond
the singleton term is

```text
p_B(z_B,m)=z_B+O(a^2)                                      (13)
```

with its first thermal derivative obeying the same bound.

### Occupied-background white polymers

An isolated matching-white hole already has activity `z_W=O(a^2)` and is
the retained second gas.  A connected or nested non-singleton defect has
at least two additional boundary units after the component reward and
contributes `O(a^3)` to the local pressure.  Interactions between retained
holes also begin at local order `z_W^2=O(a^4)`.

Thus, in a volume N,

```text
log Z_cloud
 =N[p_B+p_W]+boundary terms,
N O(a^2)=N/m^4=beta^2.                                    (14)
```

The key point is locality.  The connected correction is `N a^2`, not
`N^2 a^2`: disconnected placements exponentiate into the pressure and
must not be counted as one new polymer.

For the same-area axis/tilted pair, every contractible animal below the
injectivity radius has exactly the same embedding count per site.  Hence
the whole volume term in (14), including its thermal derivative, is a
common pressure.  It changes the pooled thermal chart but not original U,
which is invariant under a common analytic clock reparametrization.
The first geometry-dependent bulk animal must wrap the quotient and has
weight

```text
O(exp[-kappa ell log m]).                                  (15)
```

where ell is the shortest nonzero Manhattan period.

Only polymers meeting the carrier collar escape this common factor.
Their number is `O(L+S)`.  After subtracting the already retained
singleton collar, their endpoint contribution is bounded by

```text
O[(L+c)a^2].                                                (16)
```

The straight `La^2` part is common to the two rigid endpoint conventions
or may be retained as a finite collar coefficient.  The rough excess is

```text
ca^2=beta/m^3=o(beta/m).                                   (17)
```

The matching-white term has this scale already at singleton order.  All
larger connected polymers are lower.

## 4. Original thermal mark

Differentiate the polymer/collar representation with respect to the
common thermal coordinate `log h`.  A local black singleton contributes
one unit of K, a white hole minus one, and a connected animal at most its
bounded size.  Absolute convergence of the local animal series therefore
allows termwise differentiation and preserves (11), (13), and (16).

The extensive derivative of `N[p_B+p_W]` is common to the same-area
geometries.  It is absorbed into the pooled root and denominator exactly
as a common clock change.  What remains in the angular numerator is only
the differentiated carrier collar and wrapping term.  Since q and E are
bounded sector observables and the two-phase root has a nonzero thermal
slope, root motion introduces no additional N or endpoint factor.  Hence

```text
|Delta original-U mark|/|leading endpoint mark|
 <=C[beta/m+beta/m^3+exp(-kappa ell log m)]
 =o(1).                                                     (18)
```

Equation (18) is an upper bound, not a claim that the `beta/m` coefficient
is nonzero.  Occupation complement cancels its cloud-free part wordwise;
some occupied-corner collars, such as the shortest L-triomino pair, cancel
again at first local order.  The first coefficient not fixed by symmetry
is the black-singleton eligibility response of the complete rough collar.
Its natural scale is `beta/m`; white-hole and connected-polymer responses
start at `beta/m^3`.

## Scientific card

- **No hidden range state:** endpoint capillary range is `sqrt(c)`, but
  the physical eligibility collar is a tube with expected size `O(c)`.
- **Exact endpoint control:** `c partial_c log J=4c+O(1)` and has
  exponential moments; no crude c amplification is used.
- **Common pressure:** all contractible bulk polymers contribute the same
  local pressure to equal-area geometries; the first nonideal volume term
  is `beta^2` and is a common clock/normalization.
- **First geometry-dependent cloud scale:** black-singleton rough collar
  `beta/m`; matching-white holes and connected boundary polymers
  `beta/m^3`.
- **Original source:** the same orders hold after one thermal derivative;
  the relative endpoint correction is `o(1)` under beta->0.
- **Boundary:** the coefficient and sign of the complete rough-collar
  singleton response are not determined here.  No statement is made for
  a regime with beta bounded away from zero.
