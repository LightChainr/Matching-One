# Independent Burnside certificate for the bounded gadget census

For a terminal permutation, the induced action permutes the possible graph edges. A labeled graph is fixed exactly when its edge indicator is constant on every induced edge cycle. If the cycle lengths are `c_j`, the edge-refined fixed-graph polynomial is

```text
product_j (1 + x^c_j).
```

Burnside averaging over the full terminal symmetric group gives:

| terminals | group order | fixed-graph sum | graph orbits |
|---:|---:|---:|---:|
| 3 | 6 | 120 | 20 |
| 4 | 24 | 2,160 | 90 |

Coefficient-wise averaging also gives the exact edge-count histograms `1,2,4,6,4,2,1` and `1,2,5,11,17,18,17,11,5,2,1`. These agree with the canonical enumeration on `main`, although the Burnside derivation does not call graph canonicalization.

## Boundary

This certifies only total and edge-refined orbit counts for three/four terminals plus one fixed internal vertex. It does not address connectivity filters, probability polynomials, planarity, tilings, self-duality, candidate ranking, critical thresholds, or rigorous bounds. Issue #13 remains open.
