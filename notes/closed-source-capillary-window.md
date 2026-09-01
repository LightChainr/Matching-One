# The mesoscopic particle count is also the winding capillary parameter

**New obstruction and new scaling variable.**  The mesoscopic two-cloud
theorem controls the bulk rank-zero/rank-two balance, but it does not make
the first fixed-size winding shell uniform in the growing axis length.
There is an exact family of essential boundary fluctuations whose sum is

```text
F_L(m)=sum_j binomial(L,2j)m^(-2j)
      =[(1+m^-1)^L+(1-m^-1)^L]/2.                             (1)
```

Consequently the capillary variable is

```text
chi=L/m,       chi^2=N/m^2=lambda.                            (2)
```

The parameter that counts dilute rank-zero particles is therefore the
square of the parameter that roughens an axis winding interface.  A
diverging mesoscopic cloud (`lambda -> infinity`) automatically lies
outside the rigid-interface regime, even when `N/m^4 ->0` keeps the full
phase pure and makes the bulk root expansion elementary.

This separates two earlier statements that are both correct in their own
limit order:

- the two-cloud root and denominator are controlled under
  `N/m^4 ->0`;
- the fixed-L strong-source coefficient of original U is not automatically
  a uniform growing-L amplitude in that whole window.

No stochastic sample or fitted interface tension enters this conclusion.

## 1. An exact two-row family of essential curves

Take the axis `L x L` torus and one straight horizontal cut-dual boundary.
Allow it to run on either of two adjacent dual rows.  At a selected column
it changes rows through one vertical dual edge.  To return to its initial
row after one horizontal period, the number of changes must be even.

Choose any `2j` of the L column boundaries.  Between successive choices
the path follows one of the two horizontal rows, and at every choice it
switches rows.  The result is an embedded essential curve:

- it uses every horizontal column once;
- its vertical changes alternate between the two rows;
- it has no repeated edge or self-intersection;
- its length is exactly `L+2j`.

Conversely each curve in this declared two-row family is specified by its
even transition set.  Relative to the straight curve, its source boundary
weight is `m^-2j`.  Summing all transition sets proves (1).  The empty set
is the rigid straight boundary.

For a rank-one stripe of width at least three, keep its second boundary
straight and place the fluctuating boundary in the adjacent two-row slab
away from it.  This produces an actual occupied configuration in the same
rank-one homology class with

```text
g=(2L-1)+2j.                                                   (3)
```

Thus (1) is not a formal walk count detached from the lattice occupation
model.  It is a positive subfamily of the exact higher winding shells.

At the two-cloud root, different transition sets change the enclosed area
by at most L.  Their additional factor is bounded by
`exp[O(L log h)]`.  Under the pure-full gate,

```text
L log h=O(L/m^2)=O(sqrt(N/m^4)) ->0,                           (4)
```

so the area weights do not alter the asymptotics of (1).

## 2. Three capillary regimes inside the same bulk theorem

Elementary expansion of (1) gives

```text
F_L=1+binomial(L,2)m^-2+O(L^4/m^4).                           (5)
```

Hence the first omitted fixed-shell layer is already of relative order
`L^2/m^2=lambda`, not merely `m^-2` with a size-independent coefficient.
More generally, if `L,m -> infinity`,

```text
chi ->0:       F_L ->1,
chi ->c:       F_L ->cosh(c),
chi ->infinity,
L/m^2 ->0:     F_L=(1/2)exp[chi+o(chi)].                       (6)
```

The first line is the rigid-interface regime.  It is equivalent to
`N/m^2 ->0`, the zero-particle corner of the earlier Poisson theorem.  A
finite Poisson cloud gives a finite capillary dressing `cosh(sqrt(lambda))`.
The genuinely mesoscopic cloud gives exponentially many allowed two-row
fluctuations on the `exp[L/m]` scale.

This growth is still subleading to the interface tension:

```text
log F_L /(L log m) <= 1/(m log m)+o(1) ->0.                    (7)
```

So it does not erase the raw winding exponential rate.  It does erase the
claim that the coefficient of one or two fixed shells must approach its
fixed-L value throughout the larger joint window.

## 3. Consequence for the existing fixed-L U law

The completed fixed-size theorem gives

```text
Ustar/A_N
 =-(L^2-6L+6)/Delta * m^(-(2L+1))
  +O_L(m^(-(2L+3))).                                         (8)
```

Its two-power delay from the raw `2L-1` barrier is exact: reciprocal
stripe centering kills the minimum-shell thermal derivative, and root
motion, the next rank-one layer and ordinary normalization combine at
`2L+1`.

Equation (1) identifies why the subscript L on the remainder matters.
When `chi` is not small, infinitely many `g=2L-1+2j` shells have a
nonvanishing collective weight.  Their contribution must be resummed
before differentiating at the moving root.  The bulk two-cloud result
does not perform that interface resummation.

In particular, the fixed-shell normalized prediction

```text
-Delta Ustar m^(2L+1)/[A_N(L^2-6L+6)] ->1                     (9)
```

has an elementary sufficient rigid gate `L/m ->0` plus uniform control of
the remaining interface families.  It is not a consequence of
`N/m^4 ->0` alone.  In the mesoscopic regime, failure of (9) can be a
capillary dressing rather than a sign reversal, denominator collapse or
new bulk field.

The exact positive subfamily (1) by itself does not determine the sign of
the root-differentiated U after all interface families are resummed.
Accordingly this note does not upgrade the fixed-L negative sign to the
whole mesoscopic window.  It does retain the raw logarithmic tension rate
and supplies the missing variable on which the amplitude can depend.

## 4. Revised strong-source phase diagram

There are now three distinct scales:

```text
bulk particle count:       lambda=N/m^2=(L/m)^2,
full-hole count:           tau=N/m^4=(L/m^2)^2,
essential interface cost: L log m.                            (10)
```

They answer different questions.

| Regime | Bulk sectors | Axis interface |
|---|---|---|
| `lambda ->0` | empty/full | rigid; fixed shells can dominate |
| `lambda ->zeta in (0,infinity)` | finite black Poisson cloud | finite capillary dressing |
| `lambda ->infinity`, `tau ->0` | diverging black cloud/full | rough but high-tension interface |
| `tau` nonvanishing, local `N/m^6 ->0` | two particle/hole clouds | requires separate relative contour control |

The third row is the new mesoscopic coexistence theorem.  Its root remains
sharp while its shape-sensitive rank-one amplitude acquires an additional
capillary problem.  This is a mechanism separation, not a contradiction.

## 5. Next analytic target

The correct continuation is a capillary-resummed rank-one transfer, not
another finite local source or another fixed shell.  A useful target must:

1. retain both reciprocal stripe boundaries and their width zero mode;
2. sum the two-row family (1) and the other bounded-height excursions;
3. include the one-island/one-hole asymmetry responsible for the first
   nonzero thermal derivative;
4. differentiate only after inserting the two-cloud moving root;
5. compare the resulting axis term with the tilted barrier.

The last comparison is now sharp at the leading winding level.  The
[tilted shortest-class calculation](closed-source-tilted-shortest-winding-class.md)
finds two primitive classes of length `7k`, thin-cycle degeneracy
`(50k/7) binomial(7k,3k)`, and total minimum-class entropy at most its
square.  Its additional `m^-4k` cost beats that entropy for every
`m -> infinity` path considered here.  Thus the unresolved capillary
scaling function is an axis-interface problem; the tilted minimum class
cannot overturn it at leading exponential order.

The natural output is a scaling function

```text
Ustar/A_N
 =-m^(-(2L+1)) L^2 Phi(L/m)[1+o(1)],                          (11)
```

with `Phi(0)=1/Delta`.  The
[signed two-cloud transfer](closed-source-axis-signed-interface-transfer.md)
now closes this object for bounded `c=L/m`.  The exterior black-singleton
and interior white-hole activities have the same base at the exact root,
so the apparent width tilt and bridge-area response cancel.  Occupied-corner
resolution leaves two strict noncrossing bridges and gives

```text
Phi(c)=[I0(2c)^2-I1(2c)^2]/Delta.                            (12)
```

The companion [complete transfer](closed-source-axis-capillary-complete-transfer.md)
also gives the unmarked and area-refined laws

```text
G(c)=I0(2c)^2,
G(c,u)=I0(2c)^2 sinh(u/2)/(u/2).
```

It is strictly positive for every finite `c>=0`.  Thus capillary roughness
renormalizes the magnitude but cannot reverse the strong-source sign.  The
same signed complement pairing sharpens the rigid sufficient gate from
`L^2/m ->0` to `L/m ->0`, closing the formerly open `L << m <= L^2`
window.  The detailed [correction history](closed-source-axis-capillary-transfer-scaling.md)
shows why a black-only resummation spuriously creates a finite zero.

## Dependency card

- **Exact input:** the fixed closed-source boundary cost and the axis
  stripe classification.
- **New exact object:** the two-row essential-curve generating function
  (1).
- **Changed mechanism space:** bulk coexistence and winding-amplitude
  uniformity are separate; `lambda` links them as a particle/capillary
  variable.
- **No new evidence block:** this is a constructive subfamily of the same
  analytic model.
- **Closed here:** the entire rigid window `L/m ->0` and the bounded-`L/m`
  signed axis transfer.
- **Open:** unbounded `L/m`, fixed-m thermodynamics and higher tilted
  classes.  The minimum tilted class is already exponentially subordinate.

The companion [two-cloud synthesis](closed-source-mesoscopic-two-cloud-synthesis.md)
contains the bulk theorem and root expansion to which this interface
analysis is attached.

The separate [axis uniformity note](closed-source-axis-capillary-uniformity.md)
records the earlier absolute bound.  The signed transfer improves its
`L^2/m ->0` gate to `L/m ->0` and then supplies (12) at finite `L/m`.
