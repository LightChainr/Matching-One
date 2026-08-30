# H16 orientation-count no-go at fixed N=1105

This certificate resolves the bounded orientation-count question in parent
issue #74.

A primitive Gaussian layer with `k` distinct split-prime factors
`p = 1 (mod 4)` has `2^(k-1)` primitive orientation orbits modulo `D4`.
Consequently, the orbit capacity jumps by powers of two.

- `N=1105=5*13*17` is the first layer with four orbits.  Its four known
  orientations are all the available primitive `D4` directions.
- Resolving generic coefficients `H0,H4,H8,H12,H16` requires five independent
  evaluation rows.  A four-row system has rank at most four, so no fifth
  orientation can be added while keeping `N=1105` fixed.
- Requiring at least five orbits forces a fourth split prime.  The least layer
  is therefore `5*13*17*29=32045`, and the capacity jumps from four to eight.
  Exact rational elimination verifies that its first five listed orientations
  already have rank five on `H0` through `H16`.

## Boundary

The result is a counting obstruction and a formal rank certificate.  It does
not claim that `N=32045` is statistically or computationally desirable, does
not construct production projector weights, and does not modify the existing
four-orientation `N=1105` estimator.  Parent issue #74 remains open.
