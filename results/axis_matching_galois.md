# Axis matching-polynomial Galois certificates

Source: `scripts/certify_axis_matching_galois.py`.
Claim level: C5 for L=2 `C4`; C5 for L=3 `S9` and L=4 `S16` via Dedekind–Frobenius
plus the Jordan transposition criterion. Silent about infinite-volume `p_c`.
Axis L=5 / PR #84 is deliberately excluded.

## Pairwise gcds over Q

```text
L2_L3: 1
L2_L4: 1
L3_L4: 1
```

## Certificates

| L | degree | group | irreducible p | primitive p | transposition p |
|---:|---:|---|---:|---:|---:|
| 2 | 4 | C4 | 3 | — | — |
| 3 | 9 | S9 | 5 | 23 | 11 |
| 4 | 16 | S16 | 5 | 19 | 31 |

## What this does not establish

- algebraicity or transcendence of square-site `p_c`;
- Galois groups of diamond polynomials or of axis L=5;
- a finite-cell closed form for the infinite threshold.
