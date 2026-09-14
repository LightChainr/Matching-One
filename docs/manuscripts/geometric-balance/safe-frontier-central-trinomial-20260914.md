# Safe homology frontiers and the central trinomial state count

2026-09-14. Structural observation from the fixed-width charge-transfer oracle.

This note separates an exact machine-checked enumeration from a still-missing bijective proof. It is not needed for the Perron free-energy theorem; it explains why the safe state space is much smaller and more rigid than a generic connectivity-partition transfer space.

## 1. Safe frontier state

Process the cylinder one horizontal row at a time. A frontier state records occupied sites in the current row, their explored connectivity, and integer horizontal deck gains. A transition is rejected immediately when two routes inside one component imply inconsistent deck gain; that inconsistency is exactly nonzero horizontal winding.

`fixed_width_charge_transfer.py` constructs this state space by exhaustive BFS from the empty frontier, with no state cap or pruning heuristic.

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

## 4. Annular/type-B interpretation: conjectural bijection

The multiplicity `binom(2r,r)/2` is half the type-B Catalan / annular noncrossing count. This is precisely the combinatorics suggested by a cut annulus with `r` occupied boundary intervals: explored connectivity can join the runs on either side of the cut, while nonzero total deck winding is forbidden.

> **Conjecture (safe-frontier type-B bijection).** For a fixed cyclic occupancy mask with `r` occupied runs, safe connectivity-plus-gain states are in bijection with one orientation class of type-B noncrossing partitions on `r` signed boundary objects. Fixing the physical deck orientation removes the other class and leaves `binom(2r-1,r-1)` states.

A proof should construct both maps explicitly: read the signed annular partition from an explored universal-cover strip, then reconstruct integer deck gains and prove every resulting state is realizable by a finite safe square-site region.

## 5. Why the NN/matching state-set equality is plausible

Matching diagonals alter which safe state can be reached in one row, but before winding is created they should not change the vocabulary of annular frontier connectivities. Every safe matching diagonal is locally contractible; after adding finitely many past-row vertices it should admit an NN detour without changing deck gains.

> **Conjecture (adjacency-independent safe vocabulary).** For every `w`, NN and matching site graphs have identical reachable safe frontier states under the gain encoding above.

The machine check through `w=8` supports this. A facewise detour proof, analogous to the local replacement in digital Alexander duality, is the natural route.

## 6. Interface to the charge-transfer calculation

The formula supplies three cheap controls. A future implementation returning a different state count has a frontier/gain convention mismatch before any eigenvalue is trusted. The common G4/G8 vocabulary means the charge criterion compares two positive kernels on one canonical state space. Finally, central-trinomial growth gives an honest complexity estimate for the transparent reference implementation.

## 7. Claim boundary

The finite counts, run-number refinement, and G4/G8 state-set equality through the checked widths are machine-checked facts. The all-width type-B bijection and all-width adjacency-independent vocabulary remain conjectures pending constructive proof.
