# A real 8x8 configuration realizes the two-insertion pair pole

**Result.** The stated four paths are exactly four disjoint occupied
NN trees. They realize a double contraction of the D4-averaged local
pair tensor with itself, while each single insertion closes to zero.
Their fixed-configuration coefficient is

\[
 \boxed{v^{18}Q^{70}\,
 \frac{Q(Q-3)(3Q^2-9Q+8)}{8(Q-1)(Q-2)}.}
 \tag{1}
\]

Its Q=1 residue is `v^18/2`; its relative coefficient at Q=4 is
`5/768`. This is a microscopic, occupation-coefficientwise obstruction
to extending the unrenormalized one-insertion tensor as an analytic
two-insertion family. It is not a proof that a homogeneous summed
partition, a globally connected correlator, or normalized global U
diverges.

Base: `7681eedd`. The insertion and its original-rank convention are
fixed in [local four-port pair insertion](local-four-port-pair-insertion.md).
No configurations were enumerated or sampled for this construction.

## 1. Coordinates, edges and hypergraph clusters

Work on `(Z/8Z)^2`, with vacant marked sites x=(0,0) and y=(3,0).
The occupied set is the union of precisely these four paths:

\[
\begin{aligned}
C_1&=\{(1,0),(2,0)\},\\
C_2&=\{(4,0),(5,0),(6,0),(7,0)\},\\
C_3&=\{(0,1),(3,1)\}\cup\{(j,2):0\le j\le3\},\\
C_4&=\{(0,7),(3,7)\}\cup\{(j,6):0\le j\le3\}.
\end{aligned}
\tag{2}
\]

Every other site is vacant. Along row zero, x and y separate C1 from
C2, including across the periodic seam between columns 7 and 0.
The only occupied sites in rows 1 and 7 are at columns 0 and 3;
their contacts with row-zero path endpoints are diagonal, not NN.
Rows 2 and 6 are separated from row zero by those intervening rows,
and the upper and lower paths have no NN contact, including across
the vertical seam. Hence there are no unintended NN mergers or edges.

| component | occupied vertices | internal NN edges | incident hypergraph edge-nodes |
|---|---:|---:|---:|
| C1 | 2 | 1 | 7 |
| C2 | 4 | 3 | 13 |
| C3 | 6 | 5 | 19 |
| C4 | 6 | 5 | 19 |
| total | 18 | 14 | 58 |

The last column is `4|Ci|-B_i`: an internal NN edge was counted
twice in the incident-degree sum and must be counted once. Each
component is a path, so

\[
 K=18,\quad B_{occ}=14,\quad C_B=4,\quad
 \beta_1=B_{occ}-K+C_B=0,\quad r=0.
 \tag{3}
\]

There are 128 NN edge-nodes in total. The mixed edge count is
`4K-2B_occ=44`, and

\[
 B_{vac}=128-4K+B_{occ}=70,\qquad
 c_H=C_B+B_{vac}=74.
 \tag{4}
\]

Thus exactly four hypergraph clusters have sizes 7,13,19,19;
the other 70 edge-nodes are singleton vacant-vacant clusters.
No rank is added by the insertion diagram: the physical q=-1 and
E=1 remain those of (2).

## 2. The eight physical ports implement the required reflection

Give C3,C1,C4,C2 the independent colours a,b,c,d respectively.
The actual port map is

| marked site | N | E | S | W |
|---|---|---|---|---|
| x=(0,0) | C3 / a | C1 / b | C4 / c | C2 / d |
| y=(3,0) | C3 / a | C2 / d | C4 / c | C1 / b |

These are eight distinct physical edge-nodes, two on each of the four
nontrivial hypergraph clusters. The remaining 70 clusters do not meet
either insertion and supply exactly `Q^70`.

Let `mathcal K2` be the C4 average of `i P_[Q-2,2] i^dagger` used
for the single-site experiment. The unaveraged ordered projector is
invariant under simultaneous exchange of its two input colours and
its two output colours; together with C4 averaging this gives D4
invariance. In particular,

\[
 \mathcal K_2(a,d,c,b)=\mathcal K_2(a,b,c,d).
 \tag{5}
\]

Each single insertion has the four-singleton exterior partition,
since the four occupied components are distinct. The other marked
site still has its ordinary vacant tensor 1. Summing the four colours
therefore gives zero by the single-insertion contraction theorem.
For the double insertion, however, the two sites share all four
external cluster colours. Their sum is

\[
 Q^{70}\sum_{a,b,c,d}
 \mathcal K_2(a,b,c,d)\mathcal K_2(a,d,c,b)
 =Q^{70}\|\mathcal K_2\|_F^2.
 \tag{6}
\]

There is no factor v at either marked site, because both remain
vacant. Multiplying by the 18 occupied-site activities gives (1).

## 3. The norm is a finite colour-contraction identity

For completeness the norm can be reduced without a colour-array
enumeration. The unordered-pair projector has entries

\[
 a_2=\frac{Q-3}{Q-1},\qquad
 a_1=-\frac{Q-3}{(Q-1)(Q-2)},\qquad
 a_0=\frac{2}{(Q-1)(Q-2)},
\]

when its two unordered pairs have intersection size 2,1,0. In the
C4 average, colour patterns 2+2 contribute `(3/8)(Q)_2 a2²` to
the squared norm, patterns 2+1+1 contribute `(3/4)(Q)_3 a1²`,
and four distinct colours contribute `(1/4)(Q)_4 a0²`.
Patterns with at least three equal colours contribute zero.
Consequently

\[
 \|\mathcal K_2\|_F^2
 =\frac38(Q)_2a_2^2+\frac34(Q)_3a_1^2
       +\frac14(Q)_4a_0^2
 =\frac{Q(Q-3)(3Q^2-9Q+8)}{8(Q-1)(Q-2)}.
 \tag{7}
\]

This is an identity for integer Q>=4 followed by its specified rational
diagram continuation. At Q=4 the norm is 5/3. The unmodified colour
weight is `Q^74`, so the relative double coefficient is
`(5/3)/4^4=5/768`. At Q=1,

\[
 \operatorname*{Res}_{Q=1}\|\mathcal K_2\|_F^2=\frac12,
 \qquad
 \operatorname*{Res}_{Q=1}
       [v^{18}Q^{70}\|\mathcal K_2\|_F^2]=\frac{v^{18}}2.
 \tag{8}
\]

The analytic nonzero common factor relating this hypergraph sum to
the declared closed-source partition does not remove this
fixed-configuration pole.

## 4. The precise obstruction, and what has not been concluded

Give the original sites independent activities `v_z` and the two
vacant-site replacements independent parameters `lambda_x,lambda_y`.
The coefficient of

\[
 \lambda_x\lambda_y\prod_{z\in A}v_z
\]

in the unnormalized partition is exactly `Q^70 ||mathcal K2||²`.
The original A is uniquely specified by this multivariate monomial;
no other occupation pattern can cancel its pole coefficient. Both
corresponding single-insertion coefficients vanish. Thus finiteness
of every singly closed local insertion does not imply coefficientwise
analyticity of its two-insertion tensor family at Q=1.

If all activities are set equal, other configurations with K=18
enter the same coefficient and cancellation has not been excluded.
A normalized two-point quantity additionally divides by the full
partition; its globally connected version subtracts products of full
one-point expectations. These operations involve other occupation
patterns, whose contractions were not evaluated here. Therefore
(8) does not establish a pole in the summed homogeneous partition,
the globally connected correlator, or global U. Nor does it specify
a unique renormalized or confluent local field. The established
obstruction is the concrete, fixed-geometry multivariate coefficient.
