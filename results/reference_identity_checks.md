# Exact finite-matching reference checks

These are smoke tests for `scripts/matched_torus_reference.py`, not new mathematical claims. The identity being tested is the finite matching relation of Mertens and Ziff (2016).

All checks below used `mpmath` with 80 decimal digits.

## Axis torus

Command:

```bash
python scripts/matched_torus_reference.py \
  --geometry axis --L 3 --p 0.37 --exact
```

The geometry has `N=9` sites and exhausts all `2^9=512` site configurations.

```text
matching_function_cluster_side:  -0.707367121826062508
matching_function_wrapping_side: -0.707367121826062508
difference: -2.10843958864610464486971481025e-81
```

## Diamond torus

Command:

```bash
python scripts/matched_torus_reference.py \
  --geometry diamond --L 2 --p 0.37 --exact
```

The quotient periods are `(L,L)` and `(-L,L)`, so the geometry has `N=2 L^2=8` sites and physical period length `sqrt(2) L`. The check exhausts all `2^8=256` site configurations.

```text
matching_function_cluster_side:  -0.7472109756747842
matching_function_wrapping_side: -0.7472109756747842
difference: -3.37350334183376743179154369641e-80
```

The point of the second check is implementation-level: the same displacement-potential union-find can track non-trivial homology on the `pi/4`-rotated periodic quotient. This provides a reference target for a future C++/GPU implementation of orientation tomography.
