# Axis winding has a capillary obstruction to the fixed-L coefficient

**Result.** The fixed-L first two winding shells remain exactly those
classified in `762dbaf4`; they are not reclassified here. Their formal
coefficients are positive and their powers are correct at every L.
They are not automatically a uniform growing-L asymptotic.

For `L=5k`, `N=L²`, `m=m_L -> infinity`, put `tau=m^-1`.
The restricted axis rank-one probability at the strong root has

\[
 P_{1,axis}=L(L-1)\tau^{2L-1}[1+o(1)]
 \tag{1}
\]

under `L/m ->0`. The more delicate derivative and original-U
coefficient are uniformly the fixed-L ones under the explicit
sufficient condition

\[
 \boxed{L^4/m^2\to0\quad\text{equivalently}\quad L^2/m\to0.}
 \tag{2}
\]

Namely, at the original strong root,

\[
 \boxed{\partial_hP_{1,axis}
 =\frac{L^2(L^2-6L+6)}2\tau^{2L+1}[1+o(1)]>0,}
 \tag{3}
\]

and for the `(5k,0)/(4k,3k)` pair

\[
 \boxed{\frac{U}{A_N}
 =-\frac{L^2-6L+6}{\Delta}\tau^{2L+1}[1+o(1)]}.
 \tag{4}
\]

Condition (2) is sufficient, not claimed optimal. A concrete one-row
capillary family proves that at least `L/m ->0` is necessary for the
unrenormalized fixed-shell amplitude to tend to one. If `L/m ->c>0`,
that family alone contributes the nontrivial factor `cosh(c)`.

Base: `410015f5505dc2d8ca0e9ac904f656a4adc9fe86`.
The exact classification and algebra are in
[`762dbaf4`](https://github.com/LightChainr/Matching-One/blob/762dbaf4/notes/closed-source-square-family-leading-law.md).
The growing dilute background is controlled by
[the mesoscopic black-gas theorem](closed-source-mesoscopic-black-gas.md).
No N25 computation, enumeration, fit or sampling is repeated.

## 1. Immutable fixed-L input

Let `gamma=2L-1`. The completed classification gives

\[
 Z_{1,axis}=\tau^\gamma[A(h)+\tau^2B(h)+O_L(\tau^4)],
 \]

\[
 A(h)=2L\sum_{w=1}^{L-1}h^{Lw},
 \tag{5}

\]

where the `2L` configurations at each width w are straight horizontal
or vertical stripes. Thus every minimum configuration has

\[
 K=Lw,\qquad B_{mix}=2L,\qquad C_B=1,\qquad r=1,
 \qquad g=2L-1.
 \tag{6}

The next shell `g=gamma+2` consists exactly of one-row boundary
excursions and a straight stripe plus one isolated black site. Its
closed sums are the previously proved `B_rough+B_iso`. In particular
no additional family is being hidden in the remainder used below.

At fixed L the separately normalized root satisfies
`h0=1+tau²+O_L(tau⁴)`. Combining root movement, the full B shell and
the one-particle normalization gives exactly

\[
 [\tau^{\gamma+2}]\partial_hP_{1,axis}(h_0)
 =\frac{L^2(L^2-6L+6)}2.
 \tag{7}

This is positive for every integer `L>=5`. The raw minimum stripe
derivative vanishes at h=1 by particle-hole centering; equation (7)
is consequently a next-shell statement, not just the derivative of
the positive coefficient in (5).

## 2. An exact capillary subfamily

Fix one oriented axis boundary and allow it to run in either of two
adjacent dual rows. At a selected transition column it changes rows by
one vertical edge. Periodicity requires an even number `2j` of
transitions. For a fixed initial sheet, every choice of `2j` distinct
columns gives a simple essential boundary of length `L+2j`, and every
such two-row boundary is recovered from its transition set. Therefore
its exact relative generating function is

\[
 \boxed{F_L(m)=\sum_{j\ge0}{L\choose2j}m^{-2j}
 =\frac{(1+m^{-1})^L+(1-m^{-1})^L}{2}.}
 \tag{8}

Choose the other stripe boundary straight and leave one unused row on
the excursion side. These are genuine connected rank-one occupation
configurations, not abstract walks. Their extra g is exactly `2j`.
Different transition sets give different occupied regions. Constant
orientation/translation factors occur in both the straight and rough
families and do not change (8).

Consequently

\[
 L/m\to c<\infty\quad\Longrightarrow\quad F_L(m)\to\cosh c.
 \tag{9}
\]

If `L/m` diverges, `log F_L=L/m+o(L/m)` along the same subfamily.
Hence a claim that the complete shell multiplier tends to one requires
at least `L/m->0`. Positivity prevents the capillary family from being
canceled inside the restricted rank-one numerator. This is the precise
nonuniformity missed by reading an `O_L(tau^4)` remainder at growing L.

The barrier exponent itself is more robust. Even when F_L is nontrivial,

\[
 \frac{\log F_L}{L\log m}\to0,
\]

so the first winding cost remains `2L-1` on a logarithmic m scale.
The obstruction changes the amplitude and its thermal marking, not the
topological barrier.

## 3. Bulk particle decorations are a common pressure on the chart

The capillary correction must be separated from the dilute black gas.
For a straight stripe of width w, isolated black components can occupy
the vacant annulus after deleting the two rows adjacent to the stripe.
Its available area is `M_w=L(L-w-2)`. With `x=h/m²`, its singleton
partition is

\[
 Z_{sing}(M_w)=(1+x)^{M_w}[1+O(N/m^4)]
 \tag{10}
\]

uniformly under `N/m^4->0`; the error is the iid probability of an
adjacent singleton pair. Divide by the full dilute pressure `(1+x)^N`.
The width-w factor becomes

\[
 h^{Lw}\frac{Z_{sing}(M_w)}{(1+x)^N}
 =\left(\frac{h}{1+x}\right)^{Lw}(1+x)^{-2L}
 [1+O(N/m^4)].
 \tag{11}

On the coexistence chart `h=h_c e^{s/N}` with
`h_c=(1-m^-2)^-1`, one has `h_c/(1+h_c/m²)=1`. Also
`L/m²->0` follows already from `N/m⁴->0`. Thus (11) is
`exp(sw/L+o(1))`, with the same width dependence as the original
finite-L chart and no factor `exp(N/m²)` left over.

Nonsingleton contractible black components have total surplus activity
`O(N/m⁴)` by the aggregate polymer bound. Therefore local **bulk**
decorations contribute the common dilute pressure and a vanishing
relative correction. They do not repair the capillary multiplier (8),
which is supported on the essential boundary itself.

## 4. A uniform marked-remainder bound

Here is a sufficient control, included to distinguish the theorem (3)
from the merely necessary obstruction (9). For a rank-one configuration
write

\[
 g=\gamma+2j.
\]

For `j<L/4`, cut its essential component at a lowest lift and record:

1. the two axis boundaries and their exceptional horizontal reversals
   and vertical steps;
2. the rooted outer contours of zero-image components;
3. the side choices between successive recorded steps.

The displacement equations imply at most `2j` exceptional steps in
each boundary. Choosing their positions and the roots of contractible
components gives, after the `2L²` straight stripe choices are removed,
at most `(C L²)^j` records. The record reconstructs the cut-edge set;
that set fixes the occupation up to the already retained side choice.
Holes only add recorded contour length and are included in the same
bound. This is an overcount, but an injective upper encoding.

For `j>=L/4`, root an oriented nonbacktracking contour and use the
standard `3^n` walk bound. Its additional `m^-2j` beats the exponential
walk count uniformly once m tends to infinity. Multiple essential
components cost at least another `2L-2` and fall in this tail.
Consequently, with

\[
 q_L=CL^2/m^2,
\]

the unmarked shells after the classified B layer have relative mass
`O(q_L²/(1-q_L))` plus an exponentially smaller tail. This proves (1)
when `L/m->0`.

For the h derivative, mark one occupied site in the same records.
The absolute mark is at most N. After separately retaining the complete
`j=0,1` algebra, the normalized marked remainder is bounded by

\[
 C L^4\tau^\gamma\frac{q_L^2}{1-q_L}
 +\tau^\gamma e^{-cL\log m}.
 \tag{12}

The leading quantity (7) is asymptotic to
`(L^4/2)tau^(gamma+2)`. Dividing (12) by it gives

\[
 O(L^4/m^2)+e^{-c'L\log m}.
 \tag{13}

The rank-zero/root and thermal-denominator expansions have smaller
relative errors `O(N/m²)+O(N/m⁴)` under (2). The companion geometry has
systole `7k`, so its rank-one tail has an additional order proportional
to L and is swallowed by the last term of (13). Equations (7), (12)
and (13) prove (3)-(4) under the sufficient condition (2).

The encoding bound is intentionally not advertised as sharp. In the
intermediate window

\[
 L\ll m\lesssim L^2,
\]

the capillary obstruction vanishes, but the absolute marked remainder
above is not small enough to certify the old cancellation-sensitive
coefficient. A sharper signed/interface transfer calculation may extend
(3) there. The current theorem does not assume that unproved
strengthening.

## Scientific boundary

What is exact for all L is the minimum exponent `2L-1`, the straight-
stripe coefficient in (5), the complete next shell, and the capillary
factor (8). What is uniform is:

* the leading restricted P1 amplitude when `L/m->0`;
* the positive derivative and original-U coefficient under the stronger
  sufficient condition `L²/m->0`.

When `L/m` has a nonzero limit, the fixed-L shell amplitude must be
renormalized by an interface partition; inserting the bare coefficient
from `762dbaf4` is not justified. This remains a finite-lattice strong-
source statement, not a continuum exponent or a universal capillary
theory.

The companion [capillary-window synthesis](closed-source-capillary-window.md)
places this sufficient uniformity gate beside the exact two-cloud root and
the tilted shortest-winding comparison.
