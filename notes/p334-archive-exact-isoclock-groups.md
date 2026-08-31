# P334: no exact isoclock pair among the fixed147 real prefixes

The completed archive contains **147 distinct full safe polynomials for147
prefixes**. Every exact isoclock class is a singleton: there are zero
nontrivial groups and zero equal-clock pairs. Consequently this archive does
not supply a pair of original real prefixes with the same complete unmarked
clock and different geometry, and no structural comparison or marked-site
reliability solve was launched.

## Source and equality rule

The source is the unchanged `results/p334-all147-prefix-clocks/full_clocks.json`
first committed at9cca7bc6: all147 already-solved real N425/second/k0=252/age10/
ell(12,-19) prefixes in the frozen cohort. The original case-study exclusions
and reused-twelve policy remain those of the existing cohort; this task adds
no configurations.

Each grouping key is the complete tuple of174 exact integer coefficients
`(f0,...,f173)`, including trailing zeros. There is no floating-point curve
comparison, tolerance, normalization change, or digest-only equality shortcut.
Since all rows have the same d=173, equality of these tuples is exactly equality
of the complete conditional uniform-insertion survival law.

| Readout | Result |
| --- | ---: |
| Completed original prefixes grouped | 147 |
| Exact polynomial classes | 147 |
| Singleton classes | 147 |
| Nontrivial classes | 0 |
| Equal-clock unordered pairs | 0 |
| New samples / network solves / marked-site solves | 0 / 0 / 0 |

The result is consistent with the previously stored complete-pair clock
census; this follow-up makes the explicit isoclock partition available as a
small archive artifact rather than treating near-matching curves as evidence.

## Scientific interpretation

The constructed five-node double-star versus C4-plus-isolate example at
250c5899 proves a general identifiability failure of the unmarked clock. It
must not be promoted to an observed degeneracy in this real-prefix archive:
none is present here. Conversely, finding147 distinct polynomials does not
prove that the clock identifies arbitrary finite graphs or unseen prefixes.

The exact clock equality after blocking marked-middle sites in the two
counterfactual interventions at0143632d is a different result. Those are
modified events, not two identical original polynomials among these147 rows.

`results/p334-archive-isoclock-groups/exact_groups.json` records every class by
its original counter and preserves the source pointer. The short scorer
`scripts/p334_archive_isoclock_groups.py` only groups saved integers; it does
not decode or solve any source network, scan a graph family, or expand the
cohort. The requested archive question is closed at this scope.
