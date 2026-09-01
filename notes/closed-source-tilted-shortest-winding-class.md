# The leading winding class on the tilted Gaussian quotient

## Result

Let

```text
Lambda_k = <P1,P2>,
P1=(4k,3k),   P2=(-3k,4k),   N=det(P1,P2)=25k^2.
```

The shortest nonzero NN winding length is `7k`.  There are exactly two
unoriented primitive homology classes attaining it, namely `+-e1` and
`+-e2` in the period basis.  Every shortest representative is a monotone,
vertex-simple, chordless NN cycle.  For either unoriented class the exact
number of such geometric cycles is

```text
D_k = N/(7k) binom(7k,3k)
    = (25k/7) binom(7k,3k).                         (1)
```

Thus the two classes contain

```text
D_tilt,thin = 2D_k = (50k/7) binom(7k,3k)           (2)
```

one-site-thick occupied minimizers.  Every one has

```text
K=7k,   g=14k-1.                                    (3)
```

The complete minimum-`g` rank-one class consists of digital annuli whose
two resolved cut-dual boundaries are shortest curves in one of those same
two homology classes.  Its exact enumeration is a noncrossing-pair problem,
not the square of (1).  Nevertheless the exact characterization gives the
uniform bounds

```text
2D_k <= D_tilt,min <= 4 binom(D_k,2) < 2D_k^2.       (4)
```

The lower bound is the thin class.  The upper bound allows either annular
side of every unordered pair in either homology class; pairs which cross,
share a cut-dual edge, or fail the alternating-face corner compatibility
are automatically overcounted.

For the equal-area axis quotient `<(5k,0),(0,5k)>`, the corresponding
numbers are

```text
ell_axis=5k,  g_axis=10k-1,
D_axis,thin=10k,
D_axis,min=2(5k)(5k-1)=10k(5k-1).                  (5)
```

In particular, the exact thin-class ratio at closed-source weight
`h^K m^-g` is

```text
Z_tilt,thin/Z_axis,thin
 = (5/7) binom(7k,3k) h^(2k) m^(-4k).               (6)
```

Even after allowing every minimum interface in (4), the tilted leading
winding class is exponentially subordinate to the axis class in the
mesoscopic joint limit `m->infinity`, `N/m^4->0`, with bounded chart
`h=(1-m^-2)^-1 exp(s/N)`.  More precisely,

```text
log[Z_tilt,min/Z_axis,min]
 <= -4k log m + 2k log alpha + o(k),
alpha = 7^7/(3^3 4^4).                              (7)
```

The conclusion concerns the leading `g` class.  It does not order all
higher-interface rank-one configurations at finite coupling.

## 1. Shortest period vectors

Every period has the form

```text
lambda(a,b)=aP1+bP2=k(4a-3b,3a+4b),   (a,b) in Z^2.
```

Its Euclidean norm is

```text
||lambda(a,b)||_2=5k sqrt(a^2+b^2).                  (8)
```

The four unit coefficient vectors give `+-P1,+-P2`, each of Manhattan
length `7k`.  If `a^2+b^2>=2`, then (8) is strictly larger than `7k`, so
its Manhattan length is also strictly larger than `7k`.  Hence these four
oriented periods are all the shortest ones.

Although the Cartesian coordinates of `P1` and `P2` have a common factor
`k`, the relevant primitive condition is inside `Lambda_k`: their period
coordinates are `e1` and `e2`.  They are therefore primitive torus
homology classes.  Up to orientation there are precisely two classes.

## 2. Monotone words are honest simple cycles

A length-`7k` lift from `0` to `P1` must use exactly `4k` east steps and
`3k` north steps.  There are

```text
binom(7k,3k)                                         (9)
```

such based words.  The `P2` class similarly uses `3k` west and `4k`
north steps and has the same count.

No projected word can repeat a vertex before closing.  A repeated vertex
would make the intervening lift displacement a nonzero period represented
by fewer than `7k` NN steps.  Nor can the projected cycle have a
nonconsecutive NN chord: the chord would split it into two shorter closed
walks; neither can carry a nonzero period, while their oriented sum carries
`Pj`.  Thus every word is a vertex-simple induced essential cycle.

This also removes two possible counting aliases.

* The cycle has a unique NN traversal up to starting point and reversal,
  because it is induced.
* A positive orientation, say `P1`, selects one traversal rather than both.

There are `N` choices of starting vertex, and every oriented geometric
cycle is counted at its `7k` vertices.  Division by `7k` proves (1).
Counting the two positive period classes proves (2).  The quotient in (1)
is automatically integral by this free starting-point action; no separate
binomial divisibility assumption is being made.

Occupying exactly such a cycle gives `Bocc=K=7k`, one occupied component,
and ambient rank one.  The closed-source identity then gives (3), in
agreement with the exact winding barrier `g=2 ell_1-1`.

## 3. What the other minimum interfaces can be

Equality in the winding-barrier proof forces all of the following:

1. one occupied component with ambient rank one;
2. no contractible occupied component;
3. exactly two essential boundary curves and no contractible boundary;
4. each essential boundary has cut-dual length `7k`.

The two resolved boundary curves have opposite orientations in the same
primitive class.  They bound a digital annulus, and the occupied component
is one of its two sides.  Conversely, a compatible resolved pair whose
chosen side is a union of site cells gives a minimum-`g` rank-one
configuration.

The word *resolved* matters.  At an alternating face, two cut-dual walks
may meet at the face centre but the occupied-NN/vacant-matching convention
rounds the two occupied corners separately.  Hence ordinary vertex
disjointness of two unsmoothed dual words is sufficient but not necessary.
The correct exact count is

```text
D_tilt,min = sum_(j=1)^2 A_j,                       (10)
```

where `A_j` is the number of side-decorated unordered pairs of shortest
dual cycles in class `ej` which admit the resolved, globally noncrossing
digital-annulus realization.  The side decoration records which annulus is
occupied and, at a shared unsmoothed face centre, the compatible corner
resolution.  This avoids silently assuming that complementing an
alternating face preserves the same resolved pair.  Equation (10) is an
exact finite definition of the remaining enumeration problem.  It excludes
self-intersecting curves automatically, since every individual shortest
word was already proved simple.

There are at most `binom(D_k,2)` unordered pairs per class and two choices
of annular side, proving the upper bound in (4).  Every thin occupied cycle
has its two compatible shortest boundaries, proving the lower bound.
Squaring the monotone path count without the `A_j` condition would count
crossing and locally incompatible pairs as physical interfaces.

On the axis quotient the shortest dual cycles are the `5k` straight rows
or columns.  Two distinct parallel curves and a choice of side give
`2 binom(5k,2)=(5k)(5k-1)` stripes per direction.  The two directions prove
the second line of (5).  This includes all widths `1,...,5k-1`; its thin
subset has `2(5k)=10k` configurations.

## 4. Entropy does not beat the extra barrier

Stirling's formula gives

```text
binom(7k,3k)
 = sqrt[7/(24 pi k)] alpha^k [1+O(k^-1)],
alpha=7^7/(3^3 4^4).                                (11)
```

Equation (6) therefore measures the exact entropy-barrier competition for
thin interfaces.  The full minimum-interface upper bound (4) can add at
most a second copy of this path entropy, while the axis denominator may be
lower-bounded by its `10k` thin cycles.

It remains to control the activity because a thick minimum annulus may
have `K` of order `N`.  In the bounded mesoscopic chart,

```text
N |log h| <= N/m^2 + O(N/m^4+1) = o(k),             (12)
```

where `N=25k^2` and `N/m^4->0` imply `k/m^2->0`.
Consequently every ratio of activity factors between two configurations
is at most `exp[o(k)]`.  Combining (4), (5), (11), and the exact barrier
difference `(14k-1)-(10k-1)=4k` proves (7).

Thus the tilted quotient's larger shortest winding length is not undone
by its exponentially many staircase words or by the unresolved annular
pair multiplicity.  The tilted leading rank-one sector is down by
`exp[-4k log m+O(k)]` relative to the axis sector.  This is a joint-limit
statement about the already named closed source, not a continuum exponent,
a claim about all rank-one states, or an enumeration of every digital
annulus.

## Scientific card

- **Exact geometry:** two unoriented primitive shortest classes, length
  `7k`, and the exact thin-cycle degeneracy (1)-(2).
- **Exact closed-source cost:** every thin cycle has `g=14k-1`; all other
  leading minimizers are precisely the compatible two-boundary annuli in
  (10).
- **Mechanism changed:** the tilted class carries path entropy
  `alpha^k`, but the extra `m^-4k` barrier wins throughout the mesoscopic
  strong-coupling limit, even after a worst-case pair-entropy allowance.
- **Still open:** an exact closed formula for `A_j`, the resolved
  noncrossing-pair count.  It is unnecessary for the leading asymptotic
  ordering and should not be replaced by an unconstrained path-pair count.
