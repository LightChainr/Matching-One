# Exact N25 joint-insertion moments for the canonical regular pair

`scripts/regular_pair_joint_u_exact.cpp` implements the fixed-origin
translation reduction for the **actual two-vertex colour closure**. It
does not evaluate a covariance of the one-vertex activation weights.
The source stays `Kreg=K2bar+K0bar`, with coefficient of K0 identically
one and insertion strength `epsilon/N` at each vacant original vertex.
This implementation note supplies no new scientific contract and does
not report an enumeration or U result. The root task owns the frozen
definition, provenance record, execution, and score.

## Inputs and command

The only accepted geometries are the original Gaussian quotient pair
`(a,b)=(5,0),(4,3)`, each with `N=a^2+b^2=25`.

```bash
c++ -O3 -std=c++17 scripts/regular_pair_joint_u_exact.cpp \
  -o /path/to/regular_pair_joint_u_exact

/path/to/regular_pair_joint_u_exact 5 0 /path/to/kernel.tsv /path/to/n25_axis.csv
/path/to/regular_pair_joint_u_exact 4 3 /path/to/kernel.tsv /path/to/n25_oblique.csv
```

All four arguments are positional. An existing readable output file is
not overwritten. Each invocation emits one CSV and a JSON completion
record on standard output, including geometry, configuration count,
kernel rows, largest absolute kernel entry, and elapsed time. There is
no Monte Carlo, seed, size scan, alpha fitting, or root solver.

The sparse kernel is the delivered
`analysis/regular_pair_spatial_kernel.tsv`, first read here from external
commit `a237968f1d7a82d26b46e83c58179dbba7f1a908`. The executor will supply
its pinned copy from `c29d8bce` and verify the existing SHA256:

```
36ae069d370b1d7a4398861c928afb41aa76885c8895c696b1bc0c97e9c314fd
```

The production kernel has 1,874 nonzero rows, with signed integer
values from `-36` to `43`. It is **not regenerated** by this producer.
The named TSV fields are `key` (also accepts `packed_key`) and `g16`.
Keys must be canonical eight-label restricted-growth encodings; the
loader rejects duplicate, malformed, or out-of-range rows. An absent
valid canonical key is exactly zero. Pin/hash verification remains in
the execution record, rather than being inferred from the number of
rows. A direct-address int32 lookup uses 64 MiB to avoid hashing in the
inner loop; all moment accumulators use signed int64.

Implementation references, all read at `a237968f1d7a82d26b46e83c58179dbba7f1a908`:

- `scripts/p337_regular_pair_activation_exact.cpp`: quotient, rollback
  NN/matching components, and original digital-Alexander q/E labels.
- `scripts/p337_regular_pair_spatial_sampler.cpp`: canonical eight-port
  encoding and signed sparse lookup semantics. Its restriction to
  separated pairs is **not** inherited: this producer includes adjacent
  vacant vertices and identifies their shared physical edge.
- `notes/regular-pair-spatial-kernel.md`: the exact normalized closure
  `g_pi=partial_Q[Q^(-|pi|) sum_colours Kreg_x Kreg_y D_pi]|Q1`,
  stored as `g16=16*g_pi`.

## Original occupation, ports, and topology

Fix `x=0` vacant and enumerate each of the `2^24=16,777,216` occupations
of the other vertices once. `K` counts occupied original vertices.
The quotient representatives and traversal match the old exact
producer. At every leaf the unmodified black NN and white matching
components give

```
q = black_components - white_components - (K - occupied_NN_edges
                                             + full_occupied_squares),
E = q^2,                         q in {-1,0,1}.
```

For every other vacant `y`, form the outside component partition in
the common oriented port order

```
(xN,xE,xS,xW,yN,yE,yS,yW).
```

An incident edge whose other endpoint is occupied receives the root ID
of that occupied NN component. An edge whose other endpoint is vacant
receives a **unique physical undirected edge ID**, from an ID range
disjoint from occupied roots. The 25-vertex torus has 50 such NN edges.
Every edge is assigned once in the positive N/E direction and the same
ID is installed at its reverse incidence. Thus, for example, if `y` is
east of x, `xE` and `yW` are the same isolated edge-node, even though
both endpoints are vacant. Other vacant-vacant edges remain distinct;
sharing only an unoccupied vertex does not join them.

Canonical labels are assigned by first occurrence across all eight
ports and packed as `sum_i label[i] << (3*i)`. The lookup value is
retained with its sign. No virtual colour equality joins the original
occupation DSU or changes its rank, `q`, `E`, or `K`.

## Integer output and the factor used by the scorer

For one fixed-origin occupation A let

```
G16_0(A) = sum_(y != 0, y vacant) g16_(0,y)(A).
```

The CSV contains rows `k=0,...,25` and exactly these columns:

| Column | Exact sum over occupations with origin vacant and K=k |
|---|---|
| `count` | number of configurations; `binomial(24,k)` for k<=24 |
| `sum_q`, `sum_e` | original q/E marginals on this origin-vacant subset, for exact accounting against the imported full baseline |
| `sum_G16` | `sum_A G16_0(A)` |
| `sum_G16_q` | `sum_A q(A) G16_0(A)` |
| `sum_G16_E` | `sum_A E(A) G16_0(A)` |
| `nonzero_G16_count` | configurations with nonzero **summed** G16_0 |
| `sum_G16_contact` | `sum_A sum_(y vacant, y NN of 0) g16_(0,y)(A)` |
| `sum_G16_q_contact` | the same contact sum weighted by original q |
| `sum_G16_E_contact` | the same contact sum weighted by original E |

The `k=25` row is all zero. The optional nonzero count is not a count
of nonzero pairs: signed pair terms may cancel within a configuration.
No division by `count`, `2^24`, the vacancy probability, or the number
of eligible y is performed. These are raw moments, not conditionally
normalized expectations or independent pair observations.

The contact split is fixed before execution: y is one of the origin's
four original NN vertices, with no learned descriptor or distance fit.
It is accumulated in the same loop as the total; no extra occupations
or pair evaluations are generated. The scorer obtains noncontact
moments by exact `total - contact` subtraction in each of the three
marked sums. Both subsets are translation invariant, so the same
`1/(16*N)` factor and full-population normalization apply separately.

For the full homogeneous source the second epsilon derivative is

```
S2(A) = (2/N^2) sum_(unordered distinct vacant x,y) g_xy(A)
      = (1/N^2) sum_(ordered distinct vacant x,y) g_xy(A).
```

The Bernoulli weight `p^K(1-p)^(N-K)` and all three marks `1,q,E`
are translation invariant on each quotient. Consequently, for any
`O` in `{1,q,E}`,

```
sum_all_A weight(A) O(A) S2(A)
    = (1/(16*N)) sum_(A:0 vacant) weight(A) O(A) G16_0(A).
```

This is the only fixed-origin factor: the scorer uses `1/(16*N)`,
with the **old full 2^25 baseline partition denominator**. It must not
normalize by the origin-vacant population. The fixed-origin procedure
does not require q/E to be independent of vacancy; the identity follows
by translating complete configurations while retaining their marks.

With the delivered bound `|g16|<=43`, every complete absolute moment
sum is bounded by `24*2^24*43=17,314,086,912`, safely below signed
int64 capacity. The implementation also has a compile-time int64 bound
covering the entire permitted int32 kernel storage range. Configuration
counts use uint64. There is no floating-point accumulation; only the
completion timer uses floating point.

## What this hands to the U analysis

The three raw moments supply the normalized q/E responses after the
usual subtraction of their baseline expectation times the unmarked
source moment. Since canonical single- and double-insertion closures
vanish at Q1, products of their unactivated coefficients and first-Q
background-weight changes do not add an extra two-site term at this
order. In particular, subtracting `Cov(a_x,a_y)` of the first-Q
one-site activation derivatives would mix in a different Q order.

The separate scorer must still propagate these responses through the
original common pooled root, its displacement, and the complete
thermal-slope normalization defining
`J2(N25)=partial_logQ partial_epsilon^2 U|_(Q=1,epsilon=0)`.
Its prospective secondary `J2_noncontact` is the same full U response
to the noncontact subset, addressing the global contact-only model's
prediction `J2_noncontact=0` without a later descriptor search. The
contact and noncontact pieces are correlated exact submoments of the
same enumerated population, not independent evidence blocks.
The producer supplies neither that score nor a predicted sign. It
does not identify a continuum field or promote conditional positive
Gram closures to an unconditional positive response.

Preparation status: source and implementation note only. Compilation
and the single formal enumeration are left to the root executor; no
enumeration, kernel regeneration, production sampling, or test grid
was run while writing this implementation.
