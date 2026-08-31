# Oriented homology C6 orbit certificate

The primitive-line certificate for Issue #156 identifies `v` with `-v`, so a
60-degree rotation acts as C3.  Before that quotient the same integer matrix

```text
R = [[1,-1],[1,0]]
```

satisfies `R^3=-I` and `R^6=I`.  It therefore acts as C6 on oriented primitive
vectors.

At the frozen hexagonal norm cutoff 13, the 36 oriented primitive vectors
partition into six exact six-cycles.  Formal C6 character sums give Gram
matrix `6I`, without floating roots of unity.

## Descent to unoriented lines

The quotient `v~-v` makes the central element `R^3` trivial.  A C6 character
of charge `q` evaluates on that element as `(-1)^q`, so it descends exactly
when `q` is even.  The descended C3 charge is `q/2 mod 3`:

```text
H4:  C6 charge 4 -> C3 charge 2
H8:  C6 charge 2 -> C3 charge 1
H12: C6 charge 0 -> C3 charge 0.
```

Odd-spin controls H1, H3 and H5 do not descend.  This states precisely which
orientation information was discarded by the earlier primitive-line oracle.

## Reproduction

```text
python3 scripts/oriented_homology_c6_orbits.py
python3 -m unittest tests/test_oriented_homology_c6_orbits.py -v
```

## Boundary

This is an exact orbit and character quotient only.  It supplies no continuum
homology-sector baseline, measured response, radial score, covariance model,
or square-site promotion decision.  Issue #156 remains open.
