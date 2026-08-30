# Square-bond transport parity theorem

The exact L=2 enumeration on main shows centered even/odd polynomial parity.
This Issue #42 continuation isolates the structural reason and avoids a larger
configuration census.

For the periodic square-bond registry, the dual edge crossing bond `i` is
another primal bond `P(i)`.  The oracle verifies for `L=2,...,8` that `P` is a
permutation, swaps horizontal and vertical orientations, and that `P^2` is the
common torus translation `(-1,-1)`.  Geometric dual transport is therefore
vacancy complement followed by `P`.  It is bijective, changes occupation
`k` to `B-k`, and swaps primal/dual wrapping up to a translation that preserves
the wrapping channels.

If an observable is even under that transport, its Bernstein aggregates obey
`a_k=a_(B-k)`; if it is odd they obey `a_k=-a_(B-k)`.  Since the Bernstein term
at `B-k` and `1/2+t` equals the term at `k` and `1/2-t`, the centered polynomial
is respectively even or odd.  Fraction fixtures exercise both implications at
every checked bond count without enumerating `2^B` states.

Run:

```text
python3 scripts/square_bond_transport_parity_theorem.py
python3 -m unittest tests/test_square_bond_transport_parity_theorem.py -v
```

## Boundary

The length list is a structural regression range, not a finite-size scaling
study.  The algebraic implication applies whenever the declared transport
parity premise holds, but it supplies no orientation amplitude, exponent, or
continuum-limit result.  Issue #42 remains open.
