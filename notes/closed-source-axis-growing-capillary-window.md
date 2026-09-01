# The growing-capillary window is controlled by endpoint interaction

## Result

Continue the signed two-cloud axis transfer of
[`closed-source-axis-signed-interface-transfer.md`](closed-source-axis-signed-interface-transfer.md).
Let

\[
 c={L\over m}\longrightarrow\infty,
 \qquad \alpha_L={c^3\over L}={L^2\over m^3}.       \tag{1}
\]

For the square axis quotient, if

\[
                         \boxed{\alpha_L\longrightarrow0}          \tag{2}
\]

then the bounded-`c` determinant law extends uniformly:

\[
 {U\over A_N}=-{L^2\over\Delta}m^{-(2L+1)}
 \left\{I_0(2c)^2-I_1(2c)^2\right\}\{1+o(1)\}.    \tag{3}
\]

Equivalently, (3) holds throughout `m >> L^(2/3)`.  In this window the
response remains negative for positive `Delta`; neither overhang entropy
nor local cloud-interface coupling can generate a sign change.

The exponent `2/3` is a uniform **endpoint-interaction** threshold, not
the threshold for ordinary one-boundary SOS corrections.  All local
changes that merely renormalize the common translation-invariant
one-boundary kernel cancel coherently in the two-path determinant and
are smaller.  At `alpha_L=O(1)`, the first unresolved classes are instead

1. width-one near-osculation packets, where the two resolved boundaries
   remain disjoint but their one-cell collars overlap at simultaneous
   turns;
2. a black exterior singleton or white interior hole whose exclusion
   collar meets such a turn/neck;
3. an overhang whose erased horizontal reversal lies in that same
   width-one interaction packet.

These are two-boundary, gap-sensitive objects.  A height-two vertical
stack or an isolated one-boundary overhang is not by itself an
`alpha_L` mechanism.

The note proves sufficiency of (2).  It does not assert that the three
endpoint classes have a nonzero limiting coefficient at fixed positive
`alpha`; determining that coefficient requires a finite-gap interacting
transfer, rather than another independent-bridge Bessel approximation.

## 1. The determinant is smaller than bulk transfer by `1/c`

Write

\[
 \mathcal I_d=I_d(2c),\qquad
 J_1=\mathcal I_0^2-\mathcal I_1^2.
\]

The standard large-argument expansion gives

\[
 \mathcal I_0^2={e^{4c}\over4\pi c}\{1+O(c^{-1})\},
 \qquad
 1-\left({\mathcal I_1\over\mathcal I_0}\right)^2
 ={1\over2c}\{1+O(c^{-1})\},                      \tag{4}
\]

and hence

\[
 J_1={e^{4c}\over8\pi c^2}\{1+O(c^{-1})\}.        \tag{5}
\]

Thus an error bounded only by

\[
 \mathcal I_0^2\,O(c^2/L)                         \tag{6}
\]

is a relative `O(c^3/L)=O(alpha_L)` error in the endpoint response.
This is the source of (1).  It is not visible in a fixed-`c` proof.

## 2. One-boundary column corrections are coherent

Before imposing nonintersection, a directed digital boundary may make
an arbitrary vertical run in one longitudinal column.  Its exact
one-column symbol is

\[
 \phi_\tau(z)=\sum_{q\in\mathbb Z}\tau^{|q|}z^q
 ={1-\tau^2\over(1-\tau z)(1-\tau z^{-1})}.        \tag{7}
\]

For `tau=c/L`,

\[
 L\log\phi_\tau(z)=c(z+z^{-1})
 +{c^2\over2L}(z^2+z^{-2}-2)
 +O\!\left({c^3\over L^2}(|z|^3+|z|^{-3}+1)\right). \tag{8}
\]

The leading term gives the Bessel kernel.  The second term is the exact
height-two-stack/colliding-column correction.  It is a common Toeplitz
perturbation of both entries in the determinant, not an arbitrary error
of size (6).

To first order put

\[
 H_d=\mathcal I_{d-2}+\mathcal I_{d+2}-2\mathcal I_d.
\]

Then the determinant perturbation from (8) is

\[
 \delta J_1={c^2\over L}
       (\mathcal I_0H_0-\mathcal I_1H_1)+O(c^4L^{-2}\mathcal I_0^2).
                                                               \tag{9}
\]

Bessel recurrence gives, uniformly as `c->infinity`,

\[
 {\delta J_1\over J_1}=-{4c\over L}
       +O(L^{-1}+c^2L^{-2}).                       \tag{10}
\]

Therefore vertical stacks remain negligible even when
`alpha_L=O(1)`, since then `c/L->0`.  The same conclusion holds for any
finite-range, height-translation-invariant one-boundary packet: it
changes the common Fourier symbol, and its bulk pressure and diffusion
renormalizations enter both determinant rows coherently.

A non-directed overhang can be erased at its first horizontal reversal,
leaving a directed bridge plus a rooted finite packet.  If that packet
does not see the other boundary, summing its vertical translations again
produces a common Toeplitz correction.  Hence isolated overhangs are not
the obstruction behind (2).

## 3. Bad packets are localized to the width-one endpoint

Use the occupied-corner resolution and erase every maximal finite
one-boundary packet into the common kernel of Section 2.  A remaining
record must contain a packet whose support meets both boundary collars.
Since the paths are ordered, this can occur only while their separation
is one or two lattice rows.  Root the first such packet at its leftmost
column.

For a directed Bessel bridge let `R` be the number of vertical cut
edges.  Differentiating `I_d(2c)` gives uniformly for `d=0,1`

\[
 E R=O(c),\qquad E R^2=O(c^2).                     \tag{11}
\]

Given the erased pair, the root of an interacting packet is determined
by either a pair of vertical events or one vertical event and one cloud
exclusion site.  The longitudinal collision probability is at most
`C R^2/L`.  The black cloud density is `a=m^-2`; the white-hole density
is `a^2`, and the two-cloud equality at `h=1+a` removes all bulk-area
terms.  The remaining collar contribution is bounded by the same
`C R^2/L` expression (the smaller terms `La` and `Na^2` are already
absorbed in it for `c>=1`).

Consequently the total interacting-packet contribution is bounded by

\[
 |E_{int}|\le C\,{c^2\over L}\,\mathcal I_0^2.     \tag{12}
\]

The erasure is injective after recording packet type, root column and
the two incoming heights; finite packet types have a geometric source
tail because each extra cut edge costs `1/m`.  Thus the constant in
(12) is uniform while `m->infinity`.

Combining (5) and (12),

\[
 {E_{int}\over J_1}=O\!\left({c^3\over L}\right)=O(\alpha_L).       \tag{13}
\]

Equations (10),(13), together with the two-cloud cancellation already
proved at bounded `c`, establish (3) under (2).

## 4. What enters when `alpha_L=O(1)`

The proof identifies the first missing state, rather than merely losing
a remainder estimate.  Let the ordered boundary heights immediately
before a column be `(y,y+1)`.  The following finite packets are not
products of two one-boundary kernels:

* both paths turn toward the common gap in the same or adjacent column;
  occupied-corner smoothing keeps the cut edges disjoint, but their
  forbidden collars overlap;
* one path makes a rooted overhang and the other path occupies the row
  needed to translate that packet freely;
* a candidate black singleton outside, or white hole inside, is removed
  from the product cloud because its NN collar meets one of those turns.

All three modify the separation-one transfer entry without making the
same modification to the bulk `d>>1` entry.  They therefore evade the
Toeplitz cancellation in (9).  Their natural total upper scale is (12),
which becomes comparable with `J_1` precisely at `alpha_L=O(1)`.

By contrast, a height-two run in one column changes every separation by
the common second harmonic in (8); calling it the first endpoint defect
would be incorrect.  Likewise a black/white cloud particle far from the
two collars is already included in the exact common bases and cannot
alter the determinant ratio.

The next object required at fixed `alpha` is therefore a finite-gap
matrix kernel on states `(separation, collar occupancy)`.  A scalar
renormalization of `c`, or one more independent-interface Bessel term,
cannot decide the correction.

## Scientific boundary

The theorem is a growing strong-source result for the axis quotient.
It proves the positive determinant mechanism uniformly for
`m >> L^(2/3)` and isolates the only packet types that can compete at
the boundary of that window.  It does not compute their fixed-`alpha`
coefficient, claim a new continuum exponent, or address fixed `m`.
