# P418 archive-first CRT de-gauging result

The exact CRT section and fixed phase masks reproduce, but the resulting mask-times-positive-Fourier cone is rejected in all four archived channels. No new Monte Carlo was run.

| channel | raw cone d2 | raw p | CRT-masked d2 | masked p | increment over raw |
|---|---:|---:|---:|---:|---:|
| plus r1 | 70.8780 | 0.3147 | 584.0221 | 0.003984 | 513.1440 |
| plus r2 | 75.6920 | 0.2151 | 668.6463 | 0.003984 | 592.9544 |
| minus r1 | 76.7658 | 0.1633 | 1135.9943 | 0.003984 | 1059.2285 |
| minus r2 | 68.0672 | 0.3785 | 1153.1817 | 0.003984 | 1085.1145 |

`0.003984 = 1/251` is the minimum attainable p value with the frozen 250-replicate bootstrap.

## Exact and replay gates

- `s_CRT(x)=405*s(x)` passes projection, all `101^2` homomorphism pairs, C4 covariance, deck annihilation and unique-fiber checks in both child groups.
- Both hands reproduce residual-phase counts `(21,20,20,20,20)` and all four attenuation witnesses from Issue #418.
- No mask value is zero; every mask Fourier spectrum is nonnegative within floating roundoff.
- The P406 raw distances and bootstrap p values reproduce with maximum error exactly zero.

## Interpretation

The known mask is algebraically valid and invertible, but positivity after CRT de-gauging is incompatible with the archived covariance. Thus the deterministic C4 phase mask cannot account for the raw-cone compatibility under the frozen ensemble-factorization and observable-transport assumptions.

Each design still has rank `69/101`; neither this rank nor an NNLS support is a physical state count. The stored spectral-mass and 101-residue bootstrap envelopes are conditional diagnostics of a rejected model, so there is no accepted CRT-degauged prediction envelope.

The next gate is to audit the one-anchor ensemble factorization and translation transport of the projective-leg root observable before assigning the residual to fields, Jordan structure or ordered memory.

