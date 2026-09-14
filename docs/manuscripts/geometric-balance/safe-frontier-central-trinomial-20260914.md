# Safe homology frontiers and no-zero-block type-B noncrossing partitions

2026-09-14. Structural observation from the fixed-width charge-transfer oracle, sharpened by a precise combinatorial match.

This note separates machine-checked enumeration from a still-missing constructive topology/combinatorics bijection. It is not needed for the Perron free-energy theorem; it explains why the safe state space is much smaller and more rigid than a generic connectivity-partition transfer space.

## 1. Safe frontier state

Process the cylinder one horizontal row at a time. A frontier state records occupied sites in the current row, their explored connectivity, and integer horizontal deck gains. A transition is rejected immediately when two routes inside one component imply inconsistent deck gain; that inconsistency is exactly nonzero horizontal winding.

`fixed_width_charge_transfer.py` constructs this state space by exhaustive BFS from the empty frontier, with no state cap or pruning heuristic.

Contract each occupied run in the current cyclic row to one frontier object.  If there are `r` occupied runs, the nontrivial state information is the way explored components connect these `r` boundary objects through the past half-cylinder, together with their relative deck lifts.

## 2. Exact enumeration through width eight

For BOTH NN (`G4`) and matching (`G8`) adjacency, the reachable safe-state counts are

| `w` | safe states |
|---:|---:|
| 2 | 3 |
| 3 | 7 |
| 4 | 19 |
| 5 | 51 |
| 6 | 141 |
| 7 | 393 |
| 8 | 1107 |

Moreover the two graphs have the same set of encoded safe states at every one of these widths. Their transition maps/weights differ; their reachable state vocabulary does not.

The count is exactly the central trinomial coefficient

\[
\boxed{T_w=[z^0](1+z+z^{-1})^w
=\sum_{j=0}^{\lfloor w/2\rfloor}\binom{w}{2j}\binom{2j}{j}.}
\tag{2.1}
\]

The implementation asserts (2.1). Its asymptotic is

\[
T_w\sim\frac{\sqrt3}{2\sqrt{\pi w}}3^w.                       \tag{2.2}
\]

## 3. Refinement by occupied runs

Fix a nonempty binary occupancy mask on the cyclic frontier. Let `r` be the number of occupied runs. A fully occupied row is not safe because its horizontal edges already create winding.

Exhaustive enumeration gives, for every fixed mask tested through `w=8`,

\[
\boxed{N_{\rm safe}(\text{fixed mask with }r\text{ runs})
=\binom{2r-1}{r-1}=\frac12\binom{2r}{r}.}                    \tag{3.1}
\]

Thus `r=1,2,3,4` give `1,3,10,35`. The number of cyclic binary masks with exactly `r` occupied runs is

\[
2\binom{w}{2r}.                                                \tag{3.2}
\]

Multiplying (3.1) and (3.2), then adding the empty state, gives exactly (2.1). So the central-trinomial agreement is refined by the number of cyclic occupied runs; it is not only a total-count coincidence.

## 4. The combinatorial class is now much more specific

A type-B noncrossing partition on signed objects `+/-[r]` is invariant under sign reversal and has at most one sign-invariant **zero block**.  The total type-B noncrossing count is

\[
\binom{2r}{r}.                                                \tag{4.1}
\]

More specifically, the number with **no zero block** and exactly `k` pairs of nonzero blocks is known to be

\[
\binom{r}{k}\binom{r-1}{k-1}.                               \tag{4.2}
\]

Summing over `k`, Vandermonde gives

\[
\boxed{
|NC_B^{\rm nozero}(r)|
=\sum_{k=1}^r\binom{r}{k}\binom{r-1}{k-1}
=\binom{2r-1}{r-1}.}                                         \tag{4.3}
\]

Equation (4.3) is **exactly** the safe-state multiplicity (3.1).  Thus the previous vague phrase “one half of type-B Catalan objects” can be replaced by a precise target:

\[
\boxed{\text{safe frontier with }r\text{ runs}
\quad\leftrightarrow\quad
NC_B^{\rm nozero}(r).}                                       \tag{4.4}
\]

The no-zero-block condition also has the correct topology: a sign-invariant/zero block is the natural finite signed-partition image of an explored component whose lift meets one of its deck translates, i.e. a component that already carries horizontal homology and should have been rejected from the safe state space.

## 5. Proposed explicit universal-cover map

Label the `r` occupied runs in cyclic order.  Lift the processed cylinder to an infinite strip and choose one fundamental copy of every frontier run.  For each explored component touching the frontier:

1. collect the chosen frontier runs whose selected lifts belong to one lifted component;
2. also record which frontier-run lifts in the adjacent deck copy belong to that same lifted component;
3. encode the latter as the signed partners of the former.

Deck translation sends every lifted component to another component and induces block negation.  If the explored state is safe, no lifted component is fixed by deck translation, so no block can be a zero block.  Planarity/noncrossing of disjoint explored components in the cylinder should give the type-B noncrossing condition.

This gives a natural injective candidate map

\[
\Phi:\{\text{safe frontier states on a fixed }r\text{-run mask}\}
\longrightarrow NC_B^{\rm nozero}(r).                        \tag{5.1}
\]

The integer gain stored by the transfer state is exactly the data telling whether another boundary run is met in the same fundamental lift or in the neighbouring deck lift.

## 6. What is still missing for a theorem

Two points require a constructive proof rather than count matching.

### 6.1 No hidden higher-deck ambiguity

One must prove that a safe planar frontier component cannot require additional independent data beyond the signed two-copy type-B picture.  Equivalently, after choosing an appropriate lift of each connected component, all frontier contacts relevant to the canonical state must be represented by the signed annular diagram without losing possible gain information.

The finite BFS strongly suggests this is automatic from noncrossing plus absence of a deck-invariant component, but it should be proved geometrically.

### 6.2 Surjectivity / lattice realizability

Given an arbitrary no-zero-block type-B noncrossing partition, construct a sufficiently deep explored square-lattice half-cylinder realizing it without creating horizontal homology.  The noncrossing diagram supplies disjoint topological arcs; the remaining work is to route those arcs on the discrete lattice while respecting the declared occupied frontier mask.

If both steps are proved, (4.3) immediately yields the all-width central-trinomial formula rather than merely explaining it numerically.

## 7. NN versus matching vocabulary

Matching diagonals alter which safe state can be reached in one row, but before winding is created they should not change the annular connectivity vocabulary.  Every safe matching diagonal lies in a contractible local face configuration and should admit a finite NN detour in the explored past without changing the signed annular partition/deck gains.

This leads to the all-width conjecture

\[
\boxed{\mathcal S^{safe}_{4,w}=\mathcal S^{safe}_{8,w}.}       \tag{7.1}
\]

The machine check through `w=8` supports (7.1).  Once the type-B bijection is proved separately for both adjacencies, (7.1) follows immediately because both are then canonically identified with the same `NC_B^{nozero}` object for each occupancy mask.

## 8. Interface to the charge-transfer calculation

The combinatorial identification supplies several cheap controls.

- A future implementation returning a different state count has a frontier/gain convention mismatch before any eigenvalue is trusted.
- The common G4/G8 vocabulary means the charge criterion compares two positive kernels on one canonical signed-partition state space.
- Central-trinomial growth gives an honest complexity estimate for the transparent reference implementation.
- A successful type-B recoding may permit a much more compact/direct transition generator than the current BFS state discovery.

## 9. Literature boundary

The total type-B noncrossing count `binom(2r,r)` is classical.  The no-zero-block refinement (4.2) appears explicitly in the type-B noncrossing literature and sums to `binom(2r-1,r-1)`.  What is not taken from that literature is the proposed identification with safe percolation frontier connectivities; that remains the repository's topology/combinatorics bridge to prove.

## 10. Claim boundary

The finite counts, run-number refinement, and G4/G8 state-set equality through the checked widths are machine-checked facts.  The published combinatorial count of no-zero-block type-B noncrossing partitions is exact.  The actual bijection (5.1), its surjectivity, and all-width adjacency-independent vocabulary remain author-level conjectures pending constructive proof.
