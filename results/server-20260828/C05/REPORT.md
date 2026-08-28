# C05 threshold-rank Newman--Ziff reference and CPU pilot

The reference implementation freezes the off-by-one convention as follows:

- `K_plus`: first 1-based occupation rank where the black primal graph has a
  rank-2 wrapping component;
- `K_minus`: first rank where the complementary white matching graph no longer
  has a rank-2 wrapping component;
- a reverse white sweep with first-cross rank `r` maps to `K_minus=N-r+1`.

Tiny all-permutation regressions reproduce direct Bernstein evaluations of
`M(p)` and its analytic first derivative to more than 60 decimal digits. Every
sample satisfies `K_minus <= K_plus`.

The Huawei ARM pilot used an axis `L=8` torus and 100,000 deterministic
counter-keyed permutations. It retained integer marginal histograms, the
sparse joint histogram, and integer joint moments. Results:

```text
mean K_minus:       35.59763
mean K_plus:        41.10408
mean rank gap:       5.50645
rank-gap SD:         4.21495
reconstructed root:  0.5925842499338915123
M(p_ref):             0.00134850171969969
M'(p_ref):            8.33469658658750
wall time:           83.85 s
```

The compact histogram contract is therefore validated. The pure-Python
reference is not an efficient million-sample, multi-size production engine:
this pilot projects to roughly 14 minutes per million samples at `L=8` on one
server process, with steeper growth at larger sizes. Production should port
the frozen convention to C++/GPU or use explicitly disjoint counter chunks;
blindly scaling the Python implementation is not recommended.
