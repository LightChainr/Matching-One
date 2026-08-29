# Exact tiny controls for the open-boundary certification event

Status: semantics oracle for Issue 112. No production probability or threshold bound is claimed.

## Frozen geometry

The domain is an open `2s x s` rectangle with no periodic identification. Its left and right
`s x s` halves are the two cells of one renormalized bond. Vertices are ordered row-major with
`y` outermost and `x` innermost.

Two microscopic graphs are supported:

- `square`: horizontal and vertical nearest-neighbor site adjacency;
- `matching`: the square edges plus both diagonals of every unit face.

The second graph is the square site-matching lattice needed for the lower side of a confidence interval.

## Frozen event

For a given open-site configuration:

1. compute connected components separately inside each half;
2. reject an empty half;
3. reject a half if two or more components tie for maximum site count;
4. select the unique largest component in each half;
5. recompute components in the full rectangle and accept exactly when the two selected components
   belong to one union component.

Selection occurs inside each half before union connectivity is tested. This ordering is part of the
event definition and is locked by tests.

## Exact enumeration

For `n=2s^2` sites, define coefficients `c_k` by

```text
P_p(E) = sum_k c_k p^k (1-p)^(n-k).
```

All `2^n` configurations give:

| graph | s | edges | coefficient vector `c_0,...,c_n` | successes/configurations | `P_1/2(E)` |
|---|---:|---:|---|---:|---:|
| square | 1 | 1 | `[0,0,1]` | 1/4 | `1/4` |
| matching | 1 | 1 | `[0,0,1]` | 1/4 | `1/4` |
| square | 2 | 10 | `[0,0,2,8,19,24,20,8,1]` | 82/256 | `41/128` |
| matching | 2 | 16 | `[0,0,4,20,41,44,26,8,1]` | 144/256 | `9/16` |

In the first two rows “successes” is one successful configuration out of four; the probability column
is the same fraction. At `s=2`, the matching graph has no largest-cluster ties inside a half because
one `2 x 2` matching cell is complete on its four vertices.

## Nonmonotonicity is real and must not be hidden

The square NN unique-largest event is not increasing under opening sites. At `s=2` the exact oracle
finds 28 success-to-failure single-site additions.

The smallest frozen example uses row-major indices:

```text
mask 6:  open (1,0),(2,0)              -> success
open index 4=(0,1), enlarged mask 22   -> left_largest_tie -> failure.
```

The two left sites `(1,0)` and `(0,1)` are diagonal and disconnected in the square graph, creating
two size-one largest clusters. In the matching graph the diagonal joins them, so the same addition
remains successful. The `s=2` matching enumeration has zero success-to-failure additions, but no
general monotonicity theorem is claimed for larger cells.

This warning does not invalidate the fixed-parameter one-independent reduction or its binomial test:
those require locality, independence of nonincident block bonds, the path-lifting property, and a
certified event probability at the tested parameter. It does forbid silently treating this event as a
standard increasing crossing observable during design or coupling arguments.

## Boundary and next step

The oracle establishes exact event semantics only for `s=1,2`. Production still requires:

- a scalable open-boundary component algorithm implementing the same tie rule;
- tiny-oracle differential tests against that implementation;
- fresh independent final trials under the already frozen error budget;
- no reuse of exploration data for the final binomial certificate.

The parent issue remains open for those engineering and production decisions.
