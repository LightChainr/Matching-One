# A third registry separates H4 from H12

The existing axis-minus-diagonal landing mark in Issue #121 samples cosine
harmonics at angles `0` and `pi/4`.  H4 and H12 both respond as `(1,-1)`, so
those two means cannot distinguish the harmonics.

The exact difference at a candidate third angle is

```text
cos(4 theta)-cos(12 theta) = 2 sin(8 theta) sin(4 theta).
```

Thus the alias persists exactly when `theta` is an integer multiple of
`pi/8`.  In particular, adding another axis, diagonal or half-diagonal angle
does not help.

An exact separating witness is `theta=pi/12`:

```text
               H4    H12
theta=0         1      1
theta=pi/4     -1     -1
theta=pi/12   1/2     -1
```

The extended response matrix has rank two.  The minor using the first and
third rows is `-3/2`, so the separation does not rely on floating trigonometry.

## Reproduction

```text
python3 scripts/landing_registry_spin_separation.py
python3 -m unittest tests/test_landing_registry_spin_separation.py -v
```

## Boundary

This is an angular sampling design certificate.  It does not show that an
exact `pi/12` landing registry is cheap or unbiased on the finite square
lattice, run a simulation, derive an arm exponent or OPE coefficient, or
identify `x=21/4`.  Issue #121 remains open.
