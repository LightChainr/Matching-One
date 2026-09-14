# Balanced ternary / parenthesis encoding for safe site frontiers

2026-09-14.  Constructive combinatorics programme suggested by the exact match with the dilute periodic TL zero-string basis.

## 1. The target code is extremely simple

Let

\[
\mathcal W_w
=\left\{\sigma\in\{\circ,(,)\}^w:
\#(=\#)\right\}.                                             \tag{1.1}
\]

Then

\[
|\mathcal W_w|
=[z^0](1+z+z^{-1})^w=T_w,                                    \tag{1.2}
\]

the central trinomial coefficient.

Jacobsen's dilute periodic TL zero-string basis begins at width two as

```text
circle circle,  (),  )(
```

and has generating function exactly (1.2).  The periodic basis allows `)(` because an arc may cross the declared periodic cut without being a noncontractible closed loop.

The transparent site safe-frontier state count is the same `T_w`.  This suggests a direct all-width code

\[
\boxed{\mathcal S^{safe}_w\leftrightarrow\mathcal W_w.}       \tag{1.3}
\]

## 2. Why occupied-run boundaries are the right positions

For a cyclic binary occupancy mask `b_0,...,b_{w-1}`, let

\[
\partial b=\{i:b_i\ne b_{i-1}\}                              \tag{2.1}
\]

be its set of occupied/vacant transition edges.  If the mask has `r` occupied runs, then

\[
|\partial b|=2r.                                              \tag{2.2}
\]

All safe frontier connectivity beyond the mask itself lives in the explored past and can only distinguish how these `2r` run boundaries are connected through the annulus.

The machine enumeration says a FIXED mask has

\[
\frac12\binom{2r}{r}                                         \tag{2.3}
\]

safe states.

Now fix only the transition set `partial b`, not which complementary binary mask is called occupied.  Exactly two masks have that transition set: `b` and `1-b`.  Their combined number of safe states is therefore

\[
2\times\frac12\binom{2r}{r}=\binom{2r}{r}.                  \tag{2.4}
\]

But this is exactly the number of ways to place `r` left parentheses and `r` right parentheses on the `2r` transition positions, putting `circle` everywhere else.

Thus the counts match **locally for every fixed transition-edge set**, not only after summing over masks:

\[
\boxed{
\{\text{safe states for the complementary pair }b,1-b\}
\longleftrightarrow
\{\text{balanced parenthesis signs on }\partial b\}.}        \tag{2.5}
\]

Summing over all even transition sets recovers the central trinomial coefficient.

## 3. Interface-loop construction of the map

The natural map should use interfaces rather than occupied-cluster blocks directly.

Take a finite explored half-cylinder realizing a safe state.  Draw the black/white interfaces on the medial lattice.  Every transition edge in `partial b` creates one interface endpoint at the top boundary.  Since horizontal occupied homology has not yet appeared, the relevant interface connectivity is representable by contractible/through-cut annular arcs without an occupied deck-invariant cycle.

Choose the periodic cut.  Encode a top endpoint by

- `(` if its interface arc opens in the chosen reduced annular ordering;
- `)` if it closes;
- `circle` if no interface endpoint sits at that frontier position.

Noncrossing interface arcs force equal numbers of openings and closings.  An arc crossing the periodic cut produces a prefix imbalance and hence words such as `)(`; this is allowed and is exactly why ordinary Dyck words/Catalan counting is too small.

This gives the candidate

\[
\Phi_{int}:\mathcal S^{safe}_w\to\mathcal W_w.               \tag{3.1}
\]

Complementing occupied/vacant sites leaves the interface curves fixed but reverses their orientation convention, corresponding naturally to

\[
(\leftrightarrow).                                            \tag{3.2}
\]

Hence a complementary pair of masks accounts for the two global parenthesis orientations in (2.4).

## 4. Inverse construction

Given a balanced periodic parenthesis word:

1. use the standard periodic/dilute-TL rule to draw its unique reduced noncrossing annular arc connectivity;
2. choose which side of the interfaces is occupied from one reference boundary interval; global reversal gives the complementary mask;
3. thicken occupied regions on one side of the arcs;
4. route the arcs sufficiently deep in the square-lattice half-cylinder so that their discrete realizations are mutually separated;
5. read off occupied-component connectivity and deck gains at the top row.

Because the arcs are noncrossing and no noncontractible closed interface/zero-block is introduced, the resulting occupied state should be safe.

The main proof obligation is step 4: show every reduced annular arc diagram can be realized by NN SITE geometry for an arbitrary declared top mask without accidental extra contacts.  Depth is free, so a channel-routing proof should be possible.

## 5. Relation to the type-B code

Balanced periodic parenthesis words are another standard model of the annular/type-B noncrossing connectivity module.  The no-zero-block type-B partition description in `safe-frontier-central-trinomial-20260914.md` and the ternary word code here should be related by the usual conversion between signed annular blocks and reduced parenthesis arcs.

Thus the current picture is

\[
\boxed{
\mathcal S^{safe}_w
\leftrightarrow
\mathcal W_w
\leftrightarrow
NC_B^{nozero}
\leftrightarrow
\mathcal B^{dTL}_{s=0}.}                                     \tag{5.1}
\]

The ternary word is probably the most practical computational representation; the type-B partition is the cleanest block-connectivity representation.

## 6. Computational payoff if the bijection is proved

The current reference implementation discovers states by BFS and hashes connectivity/gain tuples.  A direct ternary code would remove state discovery entirely:

1. enumerate all balanced ternary words of length `w`;
2. index them by combinatorial rank;
3. implement a row update directly on the interface word;
4. reject creation of a forbidden noncontractible loop.

This should give a smaller-memory transparent transfer and may make widths beyond the current `w=8--9` controls much easier.

It would also put the Bernoulli SITE kernel directly into the conventional dilute-TL connectivity language, simplifying comparison of subleading spectra and symmetry sectors.

## 7. Exact tests already available

Any proposed map must preserve all of:

- total count `T_w`;
- fixed-mask multiplicity `binom(2r-1,r-1)`;
- width-two states `circle circle,(),)(` after the declared interface convention;
- G4/G8 equality of the safe state vocabulary through the checked widths;
- exact charge-root polynomials at `w=2,3,4` after local transition weights are assigned.

These make the bijection unusually falsifiable.

## 8. Claim boundary

The balanced-word counting and dilute-TL basis dimensions are exact.  The local count (2.4) follows from the machine-checked fixed-mask multiplicity through current widths.  The explicit interface map, all-width fixed-mask multiplicity, and discrete inverse realization remain to be proved.
