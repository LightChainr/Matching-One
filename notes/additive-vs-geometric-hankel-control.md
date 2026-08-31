# Exact additive-versus-geometric Hankel control

This finite calibration isolates the coordinate warning in Issue 400. The
same power-law profile can have full additive-translation Hankel rank while
becoming finite rank after sampling on a geometric grid.

For the rational additive sequence

```text
g_n = 1/(n+1),
```

the size-`s` Hankel block is the Hilbert matrix. Exact elimination confirms
full rank for every `s=1,...,8`; its determinant agrees with

```text
1 / product_{j=0}^{s-1} ((2j+1) binom(2j,j)^2).
```

The sequence obeys the variable-coefficient recurrence

```text
(n+2) g_(n+1) - (n+1) g_n = 0,
```

for all 32 checked steps. The full Hankel ranks simultaneously reject any
nonzero constant-coefficient recurrence of lower order on the corresponding
finite support windows.

On the geometric grid, exact controls give:

- `(1/2)^n`: rank 1;
- `(n+1)(1/2)^n`: rank 2 after the first block;
- `(1/2)^n+(1/3)^n+(1/5)^n`: rank 3 after the second block.

Thus additive Hankel rank is not invariant under changing the sampling
coordinate. Polynomial-times-exponential logarithmic partners add finite
Jordan multiplicity on the geometric grid.

## Boundary

This certificate is a finite exact rational calibration through Hankel size
8. It does not inspect P250 production data, infer a physical field count,
estimate covariance, score a held-out prediction, or supply a continuum
interpretation. Issue 400 remains open.
