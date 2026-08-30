# Exact rank-one survival certificate

This bounded certificate uses the repository's exact integer-period homology engine on the Gaussian quotient `3+i` (`N=10`). Labels are fixed by declaring site `j` to be the coset represented by `(0,j)`.

For a current rank-one occupied set `S`, let `b1` count vacant singleton insertions that raise ambient rank to two. Let `b2` count unordered pairs whose members are individually safe but jointly raise rank to two. With `q` vacancies,

```text
s2(S) = [C(q-b1,2) - b2] / C(q,2).
```

The implementation checks this identity on every supported rank-one `N=10` state. It also independently computes each survival probability by direct subset enumeration and by composing the killed one-step kernel.

| state | occupied labels | line | b1 | b2 | s1 | s2 | future exit counts |
|---|---|---|---:|---:|---:|---:|---|
| A | 0,1,2,3,4 | (0,1) | 1 | 2 | 4/5 | 2/5 | 1:24, 2:48, 3:48 |
| B | 0,1,2,3,5 | (0,1) | 1 | 3 | 4/5 | 3/10 | 1:24, 2:60, 3:36 |

Thus `(k, rank, line, H2)` agrees for A and B, while their two-step and later future laws differ. The missing two-step coordinate is the exact cooperative pair count, not an inferred temporal mode.

## Boundary

This is a finite exact calibration result. It does not modify the frozen current-k0 pilot, authorize production sampling, test scaled horizons, learn a predictive state, or imply a continuum memory field or Q4/Jordan identification. Issue #403 remains open.
