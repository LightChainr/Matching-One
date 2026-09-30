# Projected component-pair criterion for exact safe-pair synergy

2026-09-29. Author algorithm and proof for occupied-site torus homology.
This supplies the geometric enumeration behind the
[completion-pair identities](completion-pair-synergy-20260929.md).
Implementation and benchmark results belong to the parent's separate note;
this note makes no code-performance, empirical-validation, or novelty claim.
The later [complete same-prefix readout](exact-completion-pairs-readout-20260930.md)
reports the actual L512 refinement; the original probe result quoted below
is retained as the motivation, not the latest precision.

## 1. Graph, ambient rank, and transverse potential

Let V=(Z/LZ)^2, N=L^2, L>=3. Use the simple square nearest-neighbour
graph with steps ±(1,0), ±(0,1) (z=4), or add the two steps
±(1,1) for the triangular graph (z=6).
For an oriented edge u->v, let delta(u,v) be its actual lattice step,
including across a periodic seam; delta(v,u)=-delta(u,v).
It is not the difference of the two canonical coordinates in [0,L)^2.
L>=3 ensures that distinct physical edges do not collapse to parallel
edges or loops. The proof below uses this simple-graph convention.

For occupied A, let H(A) be the rational span of the winding vectors
of all cycles in its induced graph, and r(A)=dim H(A). A closed walk's
lifted displacement is Lh, where h in Z^2 is its winding vector.
Assume r(A)=1 and choose a primitive integer D=(dx,dy) spanning H(A),
with a fixed sign convention. Define the integer linear functional

    ell(x,y)=dx*y-dy*x,        ell(u->v)=ell(delta(u,v)).          (1)

This is **global ambient rank one**, not an assumption that exactly one
occupied component wraps. Every occupied component has cycle windings
in QD; several may wrap in that direction and others may have rank zero.
At least one old cycle has nonzero winding. Under any rank-one insertion,
H(A) is contained in the new one-dimensional space, so D persists.

In each occupied component C choose a root with t(root)=0 and define
t(u) by summing ell along an occupied path from the root to u. Two such
paths differ by closed walks whose projected displacements are zero:
ell(Lh)=0 for h in QD. Hence t is path-independent and

    t(u')-t(u)=ell(u->u')       on every occupied edge.           (2)

It is a scalar potential on C, even when C wraps; a full vector-valued
lift need not be path-independent. Each component has its own additive
gauge. In an extension write its potential as T(u)=t(u)+g_C.

For any B containing A, r(B)=1 iff the edge values ell admit a scalar
potential on every component of B. Indeed, such a potential exists iff
all closed-walk sums vanish, iff H(B) is contained in ker ell=QD.
The old nonzero winding prevents rank zero. Conversely, any nonzero
projected closed-walk sum gives a winding independent of D, hence rank
two, even if that walk lies in a component separate from the old wrap.

## 2. Single-site safety and the contact sign

For a vacancy v and each occupied neighbour u in component C, form

    a_v(C;u)=t(u)-ell(v->u).                                    (3)

An extended potential must satisfy T(u)-T(v)=ell(v->u), hence
T(v)=a_v(C;u)+g_C. Thus v is immediately completing iff some C gives
two different values in (3). Inconsistency prevents a potential and
forces rank two by the criterion above. Conversely, if
each incident C gives a single value a_v(C), choose T(v) arbitrarily
and set g_C=T(v)-a_v(C). All contact edges then satisfy the potential
equations. Components not incident to v retain arbitrary gauges.

Consequently, different values from *different* components do not make
v completing: their gauges are independent. No contacts, or just one
contact in each of several components, are safe. Record one a_v(C)
per distinct incident component for every safe v; repeated contacts
within a component have already been checked for consistency.
Let I(v) be that component set, so |I(v)|<=z.

## 3. Exact two-site criterion

Take distinct safe vacancies v,w and write q=T(v)-T(w). Every shared
component C in I(v) intersection I(w) imposes

    q=q_C:=a_v(C)-a_w(C).                                      (4)

If a physical lattice edge v->w exists, it additionally imposes

    q=q_edge:=-ell(v->w).                                      (5)

In particular, the edge constraint has a minus sign. Let Q(v,w) be the
set of the values in (4), together with (5) when applicable. Then

    {v,w} is a synergy edge iff |Q(v,w)|>=2.                    (6)

Here a synergy edge means r(A+v+w)=2 although both single insertions
have rank one. Empty Q or any number of identical constraints is
consistent; a singleton constraint never causes synergy.
All comparisons are integer equalities, never congruences modulo L:
reducing modulo L would erase the very winding differences being tested.

Proof. A potential on A+v+w must satisfy every equation, so conflicting
values rule it out and force rank two. Conversely, choose q to equal
their common value, or arbitrarily when Q is empty. Set T(w)=0 and
T(v)=q. For every component incident to v choose g_C=q-a_v(C); for
one incident only to w choose g_C=-a_w(C). On a shared component these
choices agree by (4). Other components retain arbitrary gauges.
Single-site safety ensures all contacts satisfy their equations, and
(5) handles the only possible edge between v,w. This constructs the
required potential and proves rank one. It also covers isolated sites,
disconnected enlarged graphs, and multiple old wrapping components.

## 4. Discover all nonlocal edges without false candidates

Fix any total ordering of occupied component identifiers. For every
safe v and every unordered pair C<C' from I(v), emit one record

    (C,C',lambda_v,v),        lambda_v=a_v(C)-a_v(C').           (7)

Group records first by (C,C'), then by the exact integer lambda value.
For each component pair, enumerate the Cartesian products of every
two *unequal-lambda* buckets, considering each bucket pair once.
Every emitted unordered vertex pair {v,w} is a true synergy edge:

    lambda_v-lambda_w
      =[a_v(C)-a_w(C)]-[a_v(C')-a_w(C')]=q_C-q_C' != 0.         (8)

The vertices need not be physically adjacent. No subsequent rank test
is needed to filter these emissions. Equal-lambda buckets are not
crossed with themselves, even if they contain many vertices.

Separately visit each unordered physical edge with two safe endpoints.
For each shared occupied component compare q_C with -ell(v->w), using
one consistent orientation. Emit {v,w} once for this visit if any
comparison disagrees. This too emits only true synergy edges.
An adjacent safe pair with no shared component passes this check;
its physical edge alone cannot supply a conflicting constraint.

Completeness follows directly from (6). Any conflict either compares
two shared-component values, detected by (8), or a shared-component
value with the physical-edge value, detected by the adjacent scan.
There cannot be two distinct physical-edge constraints because the
lattice graph is simple. These exhaust the possibilities.

Deduplicate all emissions by the canonical unordered vertex pair.
Initialize every safe vertex's degree to zero. For each first occurrence
of {v,w}, increment e and both endpoint degrees. The result is exactly
the completion graph G_A, its edge count, and every degree, including
isolated safe sites. Duplicates are repeated true edges, not false
candidates. Local contacts discover even arbitrarily separated pairs.

## 5. Gauge invariance and a sign check

Changing t on C to t+b_C sends a_v(C) to a_v(C)+b_C. Single-component
consistency and each q_C are unchanged. A component-pair key changes
lambda_v to lambda_v+b_C-b_C', the same translation for every record
in that pair, so bucket equality and inequality are unchanged.
Root choices therefore cannot change safety, edges, or degrees.
Replacing D by -D negates the potentials and all constraint values
(up to gauges), also preserving consistency and the resulting graph.

The existing [square L4 witness](completion-pair-synergy-L4-20260929.md)
B={0,2,4,5,8,12}, with v=x+4y, checks both signs by hand. Take D=(0,1),
so ell=-x. The component C={0,4,5,8,12} has t=0 on x=0 and t(5)=-1;
the isolated component C'={2} has t(2)=0. Contacts give

    (a_3(C),a_3(C'))=(1,-1),   (a_6(C),a_6(C'))=(-2,0).

Thus lambda_3=2 and lambda_6=-2 detect the nonadjacent edge {3,6}.
For the adjacent pair 6->7, a_7(C)=1 gives q_C=-3, whereas the physical
step is (+1,0), so q_edge=+1. Their mismatch detects {6,7}.
These are algebraic checks of the cited witness, not a new census.

## 6. Output-sensitive cost and its limits

A spanning-forest traversal computes components and t in O(Nz).
Scanning vacancy contacts, checking repeated-component values, and
ordering each incidence list costs O(Nz^2). There are at most

    M=sum_(v safe) binom(|I(v)|,2) <= N*binom(z,2)               (9)

records. Grouping them costs O(Nz^2) with expected constant-time exact
hash dictionaries. Use exact integer word arithmetic and sufficient
width for potentials, differences, and counts; floating tolerances
would change the criterion. Bounds here count word operations and
use expected hash costs, not measured runtime or worst-case hashing.

For a particular true edge, only pairs from I(v) intersection I(w)
can emit it, at most binom(z,2) times. Each bucket pair visited has a
nonempty Cartesian product, so its iteration cost is charged to emitted
true edges; singleton-bucket groups incur only preprocessing cost.
The adjacent scan visits O(Nz) physical pairs; merging their ordered
incidence lists costs O(z) per pair, hence O(Nz^2) in total.
Deduplication costs expected O(1) per emission and stores each edge once.
Therefore, with e=|E(G_A)|, the total expected bounds are

    time O(Nz^2+z^2*e),        storage O(Nz^2+e).                (10)

Store buckets and stream their products into the deduplication set;
do not materialize a list of all repeated emissions. For fixed z this
is O(N+e) time and storage in the stated model. It is not a guarantee
of worst-case subquadratic work: if e is quadratic, explicit output is
quadratic. No assertion of a typical e or an implementation speedup
follows from (10). Comparison-tree dictionaries would add log factors.

## 7. Same-target conditional averaging, not recursive closure

Put m=N-|A|, nu=nu_2(A), and s=m-nu. For s>0, take a uniformly random
safe site U and define Y=nu_2(A+U)-nu. The existing successor identity
gives Y=deg_G(U), hence the exact prefix-conditional measurement is

    E[Y|A]=2e/s=2e/(m-nu).                                    (11)

This replaces the safe-probe mean by its conditional expectation for
the same prefix and target. It removes within-prefix probe randomness
but leaves variation of e/s among prefixes and birth cohorts. For q
independent safe probes, in any fixed prefix cohort the variance law is

    Var(Ybar_q)=Var(2e/s)+E[Var(Y|A)/q].                        (12)

Retain the original exact (D,nu) cells, cohort definitions, support,
weights, and prefix/batch dependence. Empirical weights depend on the
prefixes, so conditioning on the full prefix sample also preserves the
weighted estimator's expectation. Apply cell-dependent consumer factors
inside cells as in the source contract. If s=0 the conditional mean is
undefined; set neither it nor its missing cohort contribution to zero.

The [latest safe-insertion result](safe-insertion-independent-block-20260929.md)
reports 140,000 new prefixes, 220,112 probes, and 1,233 nonzero increments;
its matched increment contrast is +0.00668183 +/- 0.01208706 batch SE,
unresolved. Sparsity alone does not quantify the removable variance in
(12), and this derivation supplies no new cohort estimate. Replaying old
prefix seeds with exact evaluation reuses those configurations: it is
neither independent validation nor new prospective production.

Exact degrees also give a_j=#{v safe:deg_G(v)=j}, the current full
safe-successor increment law a_j/s. This does not close the future graph:
inserting v removes v and its graph neighbours from the safe set, while
previously noncompleting triples may create new edges among survivors.
The current degree histogram does not encode that update; even the
current graph alone is not generally a recursive closure certificate.
