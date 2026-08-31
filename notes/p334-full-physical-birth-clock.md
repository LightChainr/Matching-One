# P334: full physical birth clocks, without higher-order trigger enumeration

## Outcome

The complete monotone rank-two birth event of both real N425 checkpoints is
now exactly solved. No finite-order trigger truncation remains.

| Quantity | A: counter 43042514269 | B: counter 43042505280 |
| --- | ---: | ---: |
| **Mean true birth waiting step** | **17.73237780** | **20.77877866** |
| Median / 90% birth step | 16 / 31 | 19 / 38 |
| Largest truly rank-one-safe subset | 154 | 154 |
| Mean shortening left after quartic truncation | 0.02215221 | 0.39324888 |
| Minimal quintics, obtained without enumerating quintics | 1,141 | 7,196 |

The final true mean difference is **3.046400854498077 steps**, completing

\[
10.14475552\;\text{(pairs)}\quad\to\quad5.11083837\;(\le3)
\quad\to\quad3.41749753\;(\le4)\quad\to\quad
\boxed{3.04640085\;\text{(full physical)}}.
\]

**B remains later at every nontrivial time:** its true survival and true
conditional hazard are respectively greater and smaller at every step
`k=3,...,154`. At step 155 both hazards are one; no configuration can remain
rank one after 155 further insertions. This is a statement about the whole
physical event on these two saved checkpoints, not about a surrogate graph.

## Exact finite mapping: occupied components to transverse gains

The period matrix is `P=[[425,268],[0,1]]`, with ambient line `(12,-19)`.
Its physical vector is `P*(12,-19)=(8,-19)`. Hence the integer covector

\[
q(dx,dy)=19dx+8dy
\]

vanishes exactly on that line. Every cycle already occupied in the prefix
has q=0. Since the prefix retains its nonzero ambient line under insertion,
the resulting black graph has rank two **if and only if it contains a cycle
with nonzero q**.

For each connected occupied component choose a root and a spanning-tree
displacement `p(v)` from that root. Contract the whole occupied component.
An original oriented NN edge `u -> v` becomes an edge between their roots
(or the same vacant-site variables) with gain

\[
g(u,v)=q\bigl(p(u)+(dx,dy)-p(v)\bigr).
\]

Vacant sites have `p=0`. Occupied cycles have q=0, so the contracted gains
are independent of the spanning-tree path up to vertex gauge. The sum of
gains on every closed walk equals its original transverse displacement.
Zero-gain loops and duplicate equal-gain parallel edges can be removed.

This construction preserves the entire Boolean event **simultaneously for
all 2^173 insertion subsets**. It is not a contraction fitted separately to
pair, triple or quartic outcomes.

## From the gain graph to an ordinary two-terminal network

The actual occupied essential component has root 0 in A and root 9 in B.
Delete that one contracted root K from the *fully available* gain graph.
In each remaining connected component, integrate gains to a potential phi.
Every edge satisfies `g(u,v)=phi(v)-phi(u)`: the remainder is exactly balanced.
The full gain edges and these integer gauges are retained in
`whole_event_networks.json`, so this is a finite equality certificate, not an
assumed general annular packing theorem.

For an edge `K -> v`, define its port address `a(v)=g(K,v)-phi(v)`.
An available excursion `K -> u ... v -> K` has gain `a(u)-a(v)`.
Because every cycle avoiding K is balanced, rank-two birth is therefore
equivalent to a path in the remaining occupied/selected graph connecting
two different port addresses.

In each real checkpoint there is **exactly one component with two addresses**:

| Whole-event reduction | A | B |
| --- | ---: | ---: |
| Occupied components before contracting K | 11 | 13 |
| Contracted gain edges | 379 | 388 |
| Two port addresses | -406, 19 | -171, 254 |
| Address difference | 425 | 425 |
| Random sites in the terminal block-path core | 122 | 146 |
| Other fixed occupied-component vertices in core | 10 | 12 |
| Terminal network vertices, including two terminals | 134 | 160 |
| Treewidth upper bound of this actual network | 4 | 6 |

All other components have one address and cannot create a transverse cycle,
whatever subset is inserted. Replace the two port addresses by two distinct,
always-occupied terminals s,t. Keep the remaining ordinary NN adjacency,
with occupied-component roots fixed on and vacant sites selectable. Birth
is now **exactly s-t connectivity**. Removing branches off the block-cut-tree
s-t route is also exact; those branches cannot participate in a simple s-t
path. The 51/27 random sites outside the resulting cores contribute a free
factor `(1+z)^51` / `(1+z)^27` to safety counting.

This is the desired finite annular connection representation. Its validity
for these two instances follows from the explicitly certified gains, balanced
remainder and two port addresses; no unproved general two-port claim is used.

## Exact reliability polynomial on the small-width network

Apply the actual tree decomposition, retaining both terminals in every bag.
A state stores which bag vertices are occupied and their connectivity
partition through the processed subtree. Connected terminal states are
discarded. Child states are joined only when their bag occupations agree;
their connectivity partitions are united. Edges are introduced when their
first endpoint is forgotten. A random occupied site's z weight is applied
exactly once, when that vertex is forgotten. Thus

\[
F_X(z)=\sum_{U:\ \mathrm{rank}(B_X+U)=1} z^{|U|}
      =\sum_{k=0}^{173} f_{X,k}z^k
\]

counts **all** physically safe subsets, not only those avoiding known small
triggers. Integer Kronecker encoding uses base `2^175`; all coefficients
count subsets of 173 sites and are below `2^173`, so exact integer products
and additions introduce no coefficient carries. There is no FFT rounding.

The A/B calculations use only 97/1,391 maximum live boundary states and
1,701/19,685 join pairs, taking approximately **0.01/0.21 seconds on one local
core**. As inline bookkeeping, the new full polynomials reproduce the already
known physical coefficients through k=4; no enumeration was rerun.

Then

\[
S_X(k)=\frac{f_{X,k}}{\binom{173}{k}},\qquad
h_X(k)=1-\frac{k f_{X,k}}{(174-k)f_{X,k-1}},\qquad
E[T_X]=\sum_{k=0}^{172}S_X(k).
\]

The exact integer dominance certificate is

\[
F_B(z)-F_A(z)=z^3(1+z)^{47}Q_{104}(z),
\]

where all 105 coefficients of Q are positive. Its first coefficients are
`614,64945,3398249,117250158,3000176897`; its leading coefficient is 88.
Complete coefficients and integer hazard cross-products are saved alongside
all exact rational survival and hazard values.

The first new physical coefficients are

| k | A true safe sets | B true safe sets |
| --- | ---: | ---: |
| 5 | 1,135,541,451 | 1,142,655,849 |
| 6 | 30,737,889,399 | 31,095,018,815 |
| 7 | 704,717,140,180 | 718,064,180,957 |

Subtracting the k=5 coefficient from the preceding `<=4` polynomial directly
gives **1,141 / 7,196 minimal quintics**. Their individual site lists were
not enumerated or inferred. No `C(173,5)` loop was needed.

## Scientific consequence and scope

The earlier large pair-overlap delay is real but mostly offset by higher-order
connections: a 10.14-step surrogate gap becomes a **3.05-step exact physical
gap**, without changing its sign or producing a hazard crossing. The full
event depends on an ordinary small-width two-boundary reliability network,
not on an indefinitely expanding list of minimal hyperedges. This gives a
direct bridge from local candidate insertions to the complete collective
birth clock.

The distribution is for a fresh uniform without-replacement continuation,
conditional on each fixed real prefix: N425 second orientation, seed
20260831430425, k0=252, age=10, ell=(12,-19), H2=0, b2=14770. It is not a
population law, field identification or continuum exponent, and the two
selected checkpoints are not independent production replication.

Parent: `1614a17e10997656fdf2d5520846fff2a228a5cd`. Zero new MC, no server,
one core, no higher-order subset sweep and no full-suite rerun.

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/p334_contracted_birth_network.py
/Users/lc/python-envs/research-py311/bin/python scripts/p334_full_birth_reliability.py
```
