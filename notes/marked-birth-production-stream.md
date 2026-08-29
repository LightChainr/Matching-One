# Exact marked-birth permutation-path oracle

For every microcanonical size `k`, exhaustive permutation paths agree
coefficient-for-coefficient with exhaustive absent-site insertion sums:

```text
E[(N-k) I(next site)] = E[sum_(v absent) I_v].
```

The checked vector is `(q,q2,S_active/D_active,S_inactive/D_inactive,
S_full,D_full,chi4*S,chi4*D,q*chi4*D)` with both
real and imaginary spin-four components. All residuals in the JSON are
exactly zero. Axis `L=2` supplies nonzero direct `0->2` mass and verifies
`ell=null,D=0,S=2`. Gaussian `(2,1)` verifies that `chi4` uses the lifted
Euclidean direction `P*ell`; the `(2,1)` line gives `(-7+24i)/25`.

The production path stores Horvitz absent-site sums. A canonical Russo
score multiplies them by `N/(N-k)` and convolves with `Bin(N-1,k)`.
For the common-field product, the stored root-deleted `q_before*J_D` is
completed to full-configuration `q*J_D` by adding `p*J_D`, because
`q_after=q_before+S` and `S*J_D=J_D`.
No random-root replay and no per-step all-site scan is required.

## Frozen production schema

The sparse row is one aligned-batch histogram of
`(K1,K2,site01,site12,ell,iota01,iota12,P*ell,local-H4 marks)`.
The per-`k` path retains `q,q2`, both active-primal and inactive-matching
reverse gates, their canonical full source, `chi4(ell)S/D`, `q*chi4(ell)D`,
and local-H4 `S/D`. The full matching-function source is
`S_full=(S_active+S_inactive)/2` and
`D_full=(D_active-D_inactive)/2`.

The raw sides stay in the CSV even though exact complement pairing makes
`S_active=S_inactive` and `D_active=-D_inactive`. The angular mark uses the
primitive physical lift `P*ell`, never the period-coordinate line as an
angle. The saturation index is the gcd of raw same-line winding
coefficients before primitive reduction; rank two stores zero.

The prerevealed discovery pilot is `N=65` and its exact q2 child `N=130`
at 20,000 samples each, plus max-leverage `P50 N=145` at 10,000 samples.

## Scientific card

1. MECHANISM SPACE: separates active/inactive topology, even/odd birth sources, line polarization, and local landing geometry.
2. NOT PROVED: a finite pilot cannot identify `Q4 epsilon` or establish `x=21/4` scaling.
3. OBSERVER-SECTOR-SOURCE-GEOMETRY: `A_top` | Alexander odd/even | `chi4(ell)D/S` and local H4 | Gaussian orientation pairs.
4. DEPENDENCY GROUP: all scores and covariances reuse the same counter-coupled permutation batches.
5. UPWEIGHT OBSERVATION: complement-clean q2 sign/phase transfer of connected `q-J_D4` and `gamma_D4`.
