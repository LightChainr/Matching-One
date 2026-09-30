# Birth-history memory comes from completion geometry: an exact L3 mechanism

2026-09-29. New finite geometric computation, in Draft #838. The general
[hazard identity](completion-pivotal-hazard-kernel-20260929.md) links this
calculation to the [rank-one kernel theorem](birth-kernel-markov-characterization-20260929.md).

**Outcome.** On triangular L3, direction does not remove birth-history
memory. The one-step difference `5/112` is exactly explained by different
mixtures of configurations with three or four completion sites. One binary
rank-one phase closes the count-clock process from its uniform empty start.
On square L3, unmarked label-time rank is Markov, but revealing its rank-one
direction exposes non-Markov axial sectors. Observer and clock are part
of the scientific statement, not metadata.

These are finite N=9 results, not large-size extrapolations. They do not
identify the original-U mechanism or a continuum state dimension.

## 1. A microscopic quantity with an exact transmission law

Use the existing birth engine's occupied-graph convention: square
nearest-neighbour steps, and triangular steps with additional
`+(1,1),-(1,1)` edges, on `Z^2/(3Z)^2`. Rank is the dimension over Q
of the ambient homology image, including all occupied components.
D is its primitive unoriented line while rank is one. Inclusion preserves
D until rank two; a direct `0 -> 2` path has no rank-one D.

For a rank-one occupied set A define

    nu_2(A) = number of vacant sites v with rank(A union {v})=2.

Under a uniform permutation or iid Uniform site labels, respectively,

    next-insertion exit probability = nu_2(A)/(N-|A|),
    label-time exit intensity       = nu_2(A)/(1-t).

The latter conditions on current geometry and revealed past; each remaining
label has hazard `1/(1-t)`. Occupied count is hidden and random in label
time; substituting Nt for it is not valid. The resulting transmission is

    H_count(a,k+1)=H_count(a,k)
       - E[1{J1<=a,J2>k} nu_2(A_k)]/(N-k),
    partial_t H_label(s,t)
       = -E[1{T1<=s,T2>t} nu_2(A_t)]/(1-t).

The companion note proves these formulas and expresses the logarithmic
Markov defect as the integral of the two birth-by-cutoff cohorts' mean
completion-count difference. This is a microscopic-to-observable map,
not another descriptive covariate fit.

## 2. Exact occupied-set computation

[analyze.py](../analysis/birth-completion-geometry-20260929/analyze.py)
visits all **512 subsets per lattice**, computes lifted cycles by graph
traversal, and counts completion sites. An integer dynamic programme
counts ordered prefixes by first birth. Each rank-one set at layer k has
k! such prefixes; a completion edge contributes its prefix count times
`(N-k-1)!` to the direction-marked birth table.

This is new geometric computation, not a rerun of 9! permutations.
Collapsing directions exactly recovers the old paired-birth tables.
The old files and engine remain unchanged. Their hashes, every rank-one
mask, completion sites, prefix counts and transition rows are in
[result.json](../analysis/birth-completion-geometry-20260929/result.json).

| k | Square: nu_2 (number of configurations) | Triangular: nu_2 (number of configurations) |
|---|---|---|
| 3 | 0 (6) | 0 (9) |
| 4 | 1 (36) | 2 (54), 3 (27) |
| 5 | 2 (72) | 3 (54), 4 (27) |
| 6 | 3 (48) | 3 (9) |

There are no other rank-one layers. Square has `nu_2=k-3` throughout:
at fixed insertion count, exit cannot depend on birth history. This
explains its count-clock rank-only Markov property geometrically.
Its direction set is `(1,0),(0,1),(1,1),(1,-1)`; dropping diagonals
changes the observable. Triangular L3 has only `(1,0),(0,1),(1,1)`.

## 3. Triangular memory is a completion-mixture difference

At k=5 compare early `J1<=4` with late `J1=5`, both with `J2>5`.
Early and late ordered-prefix counts are 5184 and 4536; multiply by 4!
for full-permutation counts. Four sites remain.

| Completion type | Next-step survival | Fraction among early histories | Fraction among late histories |
|---|---|---|---|
| nu_2=3 (slow) | 1/4 | 3/4 | 4/7 |
| nu_2=4 (fast) | 0 | 1/4 | 3/7 |

Thus

    E[nu_2 | early]=13/4,    E[nu_2 | late]=24/7,
    S_early=3/16,           S_late=1/7,
    S_early-S_late=-(13/4-24/7)/4=5/112.

This reconstructs the previous `(4,5,6)` witness. History selects different
current completion geometries; this is not an intervention claiming that
changing birth age causes exit.

Direction does **not** repair the memory. All three marked birth tables
are equal; each retains the two count-time violations `(3,4,5)` and
`(4,5,6)`, with conditional covariance `1/90`. The unimodular map
`(x,y)->(x-y,x)` preserves the steps, period lattice and uniform source
and cycles the three directions. Hence `H_d=H/3` for either clock.
Direction is not the missing coordinate in this finite model.

## 4. Completion types close recursively here, not just for one step

Equal exit counts do not generally imply equal non-exit successor laws.
Those laws were computed too. Every triangular configuration with the
listed `(k,nu_2)`, across directions, has exactly this successor row;
probabilities divide the counts by `9-k`.

| k, current nu_2 | Non-exit successor counts | Exit sites |
|---|---|---|
| 3, 0 | 6 of nu_2=2 | 0 |
| 4, 2 | 2 of nu_2=3, 1 of nu_2=4 | 2 |
| 4, 3 | 2 of nu_2=3 | 3 |
| 5, 3 | 1 of nu_2=3 | 3 |
| 5, 4 | none | 4 |
| 6, 3 | none | 3 |

At fixed external count k, call the smaller available completion count
slow and the other fast. This yields four states:

    rank0, rank1_slow, rank1_fast, rank2.

Rank-one rows close **strongly** across configurations. Before first birth,
the observed past is all zero; from the uniform-permutation empty start,
rank-zero configurations are uniform within their layer. Averaging their
rows therefore gives the correct **weak** rank-zero transition. Rank two
is absorbing. Together these give a time-inhomogeneous Markov realization
from that start, not merely a fit to one-step marginals.

The saved nine 4-by-4 matrices reproduce the entire paired-birth law
exactly. Direct double jumps bypass rank-one phases; off-risk identity
rows are placeholders. This is an upper construction, not minimality,
nor strong lumpability of arbitrary rank-zero preparations. Its phases
depend on external k: it is not a four-state label-clock closure with
occupied count hidden.

## 5. Square label-time Markovness depends on the observer

[exact_marked_clock.py](../analysis/birth-completion-geometry-20260929/exact_marked_clock.py)
uses exact multinomial label thinning:

    H_d(p,q)=sum_(k<=ell) Multinomial(N;k,ell-k,N-ell;p,q-p,1-q)
                           * H_count,d(k,ell).

These are **unconditional sector masses**, normalized by 9!, not by
sector visitation. For either individual axial or diagonal direction,
respectively, on square L3:

    H_axis(p,q)=3 p^3 (1-q)^3 [(1+q)^3-p^3],
    H_diag(p,q)=3 p^6 (1-q)^3.

Their sum recovers the previous unmarked, separable Markov kernel:

    2 H_axis+2 H_diag=6 p^3 (1-q^2)^3.

The axial sector instead has determinant

    H_axis(a,c)H_axis(b,b)-H_axis(a,b)H_axis(b,c)
      = -9 a^3 b^3 (1-b)^3 (1-c)^3
          * (b^3-a^3) * [(1+c)^3-(1+b)^3] < 0

for `0<a<b<c<1`. The diagonal determinant is zero. At
`(a,b,c)=(1/3,1/2,2/3)` the axial early-minus-late survival difference
is `-1084/250187`; unmarked and diagonal differences are zero.

The marked process observes rank zero before entry, `(rank1,D)` while
rank one, then an absorbing rank-two state; it does not reveal future D
before birth. Its current axial state is insufficient for its observed
history, although coarser unmarked rank is Markov. **Weak Markovness need
not survive adding an observation.**

The geometric identity makes this explicit. In the pooled square cohort,

    E[nu_2 | T1<=p,T2>q] = 6q/(1+q),

independent of p. Within an axial cohort it is

    3 - 3(1-q)(1+q)^2 / [(1+q)^3-p^3],

which depends on p. Within a diagonal cohort it is 3. Since square
`nu_2=|A|-3`, this is hidden count selection after changing clocks,
not a contradiction of constant hazard at fixed k.

Triangular marked kernels are each one third of its unmarked kernel and
retain its nonzero label determinant. All closed forms, determinants and
completion-mean identities are checked by exact polynomial arithmetic in
[marked-clock.json](../analysis/birth-completion-geometry-20260929/marked-clock.json).
The coefficient ranks (one or two) are not positive hidden-state counts.

## 6. What changes next

The finite explanation now has a concrete chain:

    birth-history selection -> completion-geometry mixture -> exit hazard -> H.

It replaces the vague proposal to try direction or birth age as a fix.
For larger sizes, the question is whether cohort hazard differences
survive scaling and what recursively controls completion geometry.
L3 closure is not assumed to extend.

A useful new independent block would capture nu_2 and persistent D at
specified rank-one times alongside paired births. First choose one clock
and geometry, then compare early/late completion hazards within direction.
The exact transmission formula defines the consumer; adding unrelated
descriptors would not answer it. A finite-lag defect requires the hazard
over an interval, not merely at its midpoint. No larger-size run is
launched here. Existing L512 pair archives lack the geometries and cannot
be relabelled as this experiment. No causal or original-U claim follows.

## Reproduce / actual work

From the repository root, with Python 3.11 and SymPy:

```bash
python3 analysis/birth-completion-geometry-20260929/analyze.py
python3 analysis/birth-completion-geometry-20260929/exact_marked_clock.py
```

The first uses the standard library and two existing exact L3 tables
for a source-convention comparison; it took about .021 seconds on local
Python 3.11.15. The polynomial transform took about .252 seconds with
SymPy 1.14.0. Runtime fields naturally vary. The geometry script was
rerun while extending its saved phase realization, not to repeat physical
permutation acquisition. No cloud, Monte Carlo, full suite or document
wording tests ran in this step.
