# Spin aliases in the axis/diagonal landing registry

The local four-arm observable associated with Issue #121 compares two landing
registries: lattice axes at angle zero and diagonals at angle `pi/4`.  This
certificate records exactly which square-symmetric harmonics that two-point
angular sampling can and cannot identify.

## Exact response

Square symmetry permits spins `s=4k`.  A cosine harmonic sampled at the two
registries has response

```text
(cos(s*0), cos(s*pi/4)) = (1, (-1)^k).
```

Consequently the axis-minus-diagonal mark is two for `s=4 mod 8` and zero for
`s=0 mod 8`.  Through spin 32 the two complete alias classes are

```text
selected:     H4, H12, H20, H28
annihilated:  H0, H8, H16, H24, H32.
```

Both sine samples vanish for every `s=4k`, so the same two registry means do
not contain a hidden quadrature that could break these aliases.

## H4/H12 no-go

The two response columns for H4 and H12 form

```text
[[ 1,  1],
 [-1, -1]].
```

Its exact rank is one and its determinant is zero.  Any linear statistic of
only the axis and diagonal registry means therefore gives the same response
to H4 and H12 up to amplitude.  Calling the mark `H4` states its lowest
allowed harmonic, not an identification of a pure spin-four sector.

## Reproduction

```text
python3 scripts/landing_registry_spin_alias.py
python3 -m unittest tests/test_landing_registry_spin_alias.py -v
```

## Boundary

This is an observable-identifiability result only.  It runs no new landing
simulation and derives no arm exponent, OPE coefficient, or `x=21/4`
identification.  Resolving H4 from H12 requires additional angular registry
information or an independent scaling/shape argument; Issue #121 remains
open.
