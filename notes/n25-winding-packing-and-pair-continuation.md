# N25 winding packing makes the existing seam histogram sufficient at generic Q

**Result.** On an honest Gaussian square torus of area N, w disjoint
essential occupied NN components with common primitive deck slope h
satisfy

\[
 \boxed{N\geq w\bigl(\|Ph\|_1+\|Ph\|_\infty\bigr)
       \geq 2w\sqrt N\,\|h\|_2.}
 \tag{1}
\]

The first inequality includes a proved supply of **w vertex-disjoint
white matching cycles**, not an assumed white multiplicity. For each
of the N25 quotients `(5,0)` and `(4,3)`, (1) forces `w<=2`; if w=2,
the primitive deck slope is an axis, and if w=1, `|u|<=2`.
Consequently the completed `(k,g,q,bad2,n_bad3)` histogram already
separates every nonzero generic `[Q-2,2]` seam contribution. There is
no new enumeration or sampling in this note.

Base: `0dda27ba`. The finite homology-decorated closure and pair
contraction are in
[finite torus pair closure](closed-source-finite-torus-pair-closure.md);
the fixed hypergraph law is in
[hypergraph RC and twist projection](closed-source-hypergraph-rc-twist-projection.md).
The already completed histogram's field meanings are those of
`scripts/p337_s4_trace_exact.cpp`: `bad2` is a boolean indicating at
least one component with an odd first deck gain, and `n_bad3` counts
components whose first deck gain subgroup is nonzero modulo 3.

## 1. The white-cycle packing lemma

Use the usual embedded square NN graph on a torus. Assume that its
unit faces and NN/NNN neighbourhoods are not aliased; both specified
N25 quotients satisfy this. Let B be the induced occupied NN graph,
and suppose its total ambient image has rank one. Its essential
components have a common primitive slope h, up to orientation:
distinct essential simple curves of different slopes intersect on a
torus, which would join their occupied components. An essential simple
cycle in an embedded graph is primitive; rank-one saturation therefore
gives the whole image `Zh` of each such component.

**Lemma.** If there are w essential components of B, the vacant
NN+NNN matching graph contains w vertex-disjoint essential cycles of
the same primitive slope, and those cycles are disjoint from B.

**Proof.** First resolve the apparent problem of matching diagonals
crossing. Form a planar embedded white graph W as follows:

* keep every white-white NN edge;
* retain a white diagonal only in a face with exactly two white
  corners, opposite each other.

Every omitted white diagonal is in a face with at least three white
corners. It can be replaced by two white NN sides in that face. The
replacement is a homotopy in the face, so W has the same components
and the same ambient homology image as the full white matching graph.
There is at most one retained diagonal in any face. W is thus genuinely
embedded and avoids B. No intersection of two drawn white diagonals is
silently treated as a vertex.

Take a sufficiently thin regular neighbourhood of each component of
B. A connected rank-one neighbourhood is an annulus with possibly
contractible holes: it has two essential boundary circles of slope h;
its other boundary circles bound discs. One way to see this is to cut
the torus along an essential simple cycle in that component. The cut
surface is an annulus, and rank one prohibits a further connection
between its two copies that would create an independent winding.
Equivalently, the classification of subsurfaces of a torus gives genus
zero and exactly two noncontractible boundary circles. Rank-zero
neighbourhoods are planar regions with contractible boundary circles.

The w essential neighbourhoods have a cyclic order. Between successive
ones, including the last and first, there is an essential complementary
annular region, possibly with holes and contractible appendages. A
black component crossing one such annulus from one essential boundary
to the other would join the two corresponding black components. It
is therefore absent. Contractible black components can remove discs
or add holes but cannot remove its essential core. This gives w
different components of the complement of B, each containing an
essential simple curve of slope h. For w=1 this is the one annular
complement of the occupied essential neighbourhood.

It remains to put each continuous core curve on actual white sites.
Put it in general position relative to the square cellulation. Every
NN edge crossed by the curve is not in B, so has a white endpoint;
choose such an endpoint at every crossing. Consecutive endpoints lie
on one unit face. Join them by their white NN or matching diagonal
edge, replacing an omitted diagonal by the white NN two-step path
just described. These replacements stay in that face away from B and
preserve the homology of the original curve. They give a white closed
walk in the same complementary component.

Decompose that walk into simple cycles in the embedded graph W. At
least one is essential. Because it lies in the annular complementary
region it has slope h; because it is an embedded simple essential
curve it is primitive. Curves obtained in different complementary
components cannot share a white vertex. The claimed w cycles follow.
This argument uses the digital matching connectivity correctly, while
obtaining the cycles in a planar subgraph. **End proof.**

## 2. Counting sites, not contour lengths

Choose one essential occupied simple cycle from each of the w black
components. A lift of each has displacement `Ph`, so its number of
unit NN steps, equivalently its number of distinct vertices on the
cycle, is at least `||Ph||_1`. The white cycles supplied by the lemma
use NN or diagonal steps, each of infinity norm at most one; each
therefore has at least `||Ph||_infty` distinct vertices. All chosen
cycles of either colour are mutually vertex-disjoint within that
colour, and black and white vertices are disjoint. Their sites fit in
the N sites of the torus, proving the first inequality of (1).

For `z=(x,y)`, put `a=max(|x|,|y|)` and `b=min(|x|,|y|)`. Then

\[
 (\|z\|_1+\|z\|_\infty)^2-4\|z\|_2^2
 =(2a+b)^2-4(a^2+b^2)=b(4a-3b)\geq0.
\]

For a Gaussian quotient,
`P=[[a,-b],[b,a]]` and `P^T P=N I`; hence
`||Ph||_2=sqrt(N)||h||_2`. This proves the second inequality of (1).

At N=25, `w sqrt(u²+v²)<=5/2`. Since h is nonzero and primitive:

\[
 w\leq2;\qquad
 w=2\Longrightarrow h\in\{(\pm1,0),(0,\pm1)\};\qquad
 w=1\Longrightarrow |u|,|v|\leq2.
 \tag{2}
\]

These are necessary packing bounds. The proof does not assert that
every slope surviving them is realized. No existence search is needed
for their use in the already completed histogram.

## 3. The mod-6 record now determines the generic pair trace

Use the first primitive deck seam. A rank-one configuration with w
essential components of slope `(u,v)` and c0 other hypergraph
components has twisted colour character

\[
 Q^{c_0}\operatorname{Fix}(\pi^u)^w.
 \tag{3}
\]

For integer Q>=4 write `d2(Q)=Q(Q-3)/2` for the dimension of the
unordered distinct-pair irrep `[Q-2,2]`. If u=0, (3) is constant.
For w=1 and `|u|=1`, it is the point-permutation character, containing
only singlet and standard. Their pair projections vanish.

There are exactly two remaining possibilities allowed by (2):

| class | geometry within rank one | existing fields | pair colour contraction |
|---|---|---|---|
| A | w=2, `|u|=1` | `q=0,bad2=1,n_bad3=2` | `Q^c0 d2(Q)` |
| B | w=1, `|u|=2` | `q=0,bad2=0,n_bad3=1` | `Q^c0 d2(Q)` |

For A, the square of the point-permutation character is the character
of `V tensor V`; `[Q-2,2]` occurs once. For B, use the exact identity
`chi_V(pi²)=chi_Sym²V(pi)-chi_Lambda²V(pi)`. The symmetric square is
the sum of the diagonal-colour copy V and the unordered distinct-pair
carrier. The alternating square contains standard and `[Q-2,1,1]`,
but not `[Q-2,2]`. Its pair coefficient is therefore also +1.
Negative u has the same fixed-colour count. These are stable integer-Q
representation identities with rational diagram continuation.

Rank-zero and rank-two pair contractions vanish by the full numerator
selection argument; in fact their partition contractions vanish too:
their characters are constant and proportional to `Fix(pi)` respectively.
Every other rank-one histogram row also vanishes by (2)-(3).

## 4. Exact continuation without another enumeration

Let `c_H=c0+w` and let A and B now denote the two sets of histogram
configurations in the table. Including the original rank factor gives

\[
 \boxed{\Psi_{[2]}(Q,v)=d_2(Q)\left[
  \sum_{A}v^K Q^{c_H-5/2}
 +\sum_{B}v^K Q^{c_H-3/2}\right].}
 \tag{4}
\]

The old histogram retains `g=2N+1-K-Sstar`. Since
`Sstar=2c_H-r+3K-2N+1`, its rank-one rows determine

\[
 c_H=(4N+1-4K-g)/2,\qquad c_0=c_H-w.
 \tag{5}
\]

Thus a row of count n contributes exactly
`n v^K d2(Q) Q^(c0-1/2)`; no missing slope, cluster count or source
field is needed for (4). At fixed v its Q=1 value and derivative are

\[
 \left. n v^K d_2(Q)Q^{c_0-1/2}\right|_{Q=1}=-n v^K,
 \qquad
 \left.\partial_Q[n v^K d_2(Q)Q^{c_0-1/2}]\right|_{Q=1}
 =-n c_0 v^K.
 \tag{6}
\]

Equation (6) is a signed analytic sector continuation, not a positive
one-colour sector probability. To return to the closed-source path one
must still use the declared `v=y Q^(3/2)` and common partition factor
`Q^(1/2-N)`, and normalize the two geometries separately. Neither the
fixed-v derivative alone nor a projector-only derivative is the entire
thermal/root response. The direct q/E pair numerators remain identically
zero; (4) supplies the previously missing generic-Q denominator sector.
This continuation does not identify a local four-leg field or alter the
regular-endpoint identity `ell P_[2](Q)=0`.
