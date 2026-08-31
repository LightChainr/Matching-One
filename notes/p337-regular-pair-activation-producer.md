# Regular-pair completion: fixed-origin producer staged, not run

**Execution update:** the frozen producer has since completed both
once-only traversals after the theory gates and GO. See
[raw receipt](../results/p337-regular-pair-activation/README.md) and
[completed mixed response](regular-pair-interaction-result.md).
The remainder records the pre-run interface and its original status.

Base: `7681eedd`. Code: `scripts/p337_regular_pair_activation_exact.cpp`.
This file prepares the one named Q1 derivative source of the proposed
regular completion `Kreg=K2bar+K0bar=avg i(I-P1)i†`. Its actual tensor
derivation and the root's frozen contract must be accepted before GO.
Only source preparation and a compiler syntax check are performed now.

## Outside edge-node classification

The origin is representative vertex0. If it is occupied, the source is0.
Otherwise its four incident NN edge-nodes are labelled in physical
`N,S,E,W` order:

- An occupied neighbor contributes its occupied NN component root.
- A vacant neighbor contributes a unique singleton edge-node label,
  distinct from every DSU root and every other vacant-neighbor edge.

The second case is essential: it is neither discarded nor queried through
an inactive black DSU parent. With origin vacant, an edge whose other
endpoint is also vacant is an isolated equality-node in the hypergraph.
Occupied-neighbor connectivity is queried in the unchanged full outside
occupied NN graph.

The requested provisional source table, stored as integer `a4=4*a_x`, is:

| Outside partition | Arrangement | a4 |
|---|---|---:|
| 1+1+1+1 | four distinct | 4 |
| 2+1+1 | the pair is NS or EW | 4 |
| 2+1+1 | adjacent pair | 2 |
| 2+2 | NS\|EW | -2 |
| 2+2 | NE\|SW or NW\|ES | -1 |
| 3+1 or4 | any | 0 |

The implementation counts equal unordered port pairs:0,1,2,3,6 identify
these partition types. This is a four-port query only; it does not scan
all25 sites. The table is not independently declared proved by the producer.

## Raw interface

Output per geometry has K=0,...,25 and

```text
k,count,sum_q,sum_e,sum_a4,sum_qa4,sum_ea4
```

The last three sums are divided by4 exactly once by the future scorer.
All configurations are retained, including origin-occupied source zeros;
there is no population conditioning and no additional binomial factor.
By translation invariance, these first-source joint moments with K,q,E
equal the moments of the uniform site-average source. They do not give
multiple-insertion moments or a site-sum source without its explicit N
normalization.

Geometry `(5,0)/(4,3)`, component rollback, white matching/Alexander q,
NN edge counting, plaquette counting and binary traversal are inherited
from the existing producer. No seam state, new root, old Sstar/Bvac source,
support-radius scan or generalized engine has been added.

After the root's contract and GO, the fixed CLI will be `a b output.csv`.
The future two raw files and receipts belong only in
`results/p337-regular-pair-activation/`. No enumeration or score exists yet.
