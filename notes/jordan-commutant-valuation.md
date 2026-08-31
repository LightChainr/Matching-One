# Jordan commutant valuation and exact image rank

For one rank-`r` nilpotent Jordan block `N`, every endomorphism commuting with
`N` is an upper-Toeplitz polynomial

```text
A = a_0 I + a_1 N + ... + a_(r-1) N^(r-1).
```

Let `nu` be the first index with `a_nu != 0`.  Then `A=N^nu B`, where `B`
has nonzero constant coefficient and is therefore invertible in the truncated
polynomial algebra.  Acting on the top vector leaves exactly

```text
r - nu
```

nonzero Jordan-chain steps.  The zero polynomial has image rank zero.

This refines the earlier bottom-survival certificate for Issue #216:
`a_0!=0` is precisely the full-rank case `nu=0`.  It also quantifies every
possible collapse rather than only providing one negative control.

The checked controls at rank five realize image ranks 5, 4, 2 and 0.  A
brute-force test independently scans every 2x2 and 3x3 matrix with entries in
`{-1,0,1}` and verifies that every commuting matrix reconstructs from its
unique polynomial coefficients.  The existing thermal Q4 label action is the
rank-two identity polynomial and retains image rank two; its nonzero Gram norm
remains 4930.

## Reproduction

```text
python3 scripts/jordan_commutant_valuation.py
python3 -m unittest tests/test_jordan_commutant_valuation.py -v
```

## Boundary

This adds no Virasoro quotient calculation beyond the existing Q4 norm.  It
does not establish lattice overlap, fix the logarithmic coefficient, derive a
torus Ward response, or identify `P4[S']` as a top-field readout.  Issue #216
remains open.
