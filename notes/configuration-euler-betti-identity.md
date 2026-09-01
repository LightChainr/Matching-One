# Configuration Euler–Poincaré / Betti lift of the matching identity

Status: C5 finite identity for issue #111. Variance-reduction branch **closed**.

## Exact configuration identity

Represent an occupied square-torus configuration as the periodic cell complex
with occupied vertices `V`, occupied nearest-neighbour edges `E`, and occupied
elementary faces `F0`. Let `beta0_black` be the number of occupied primal
components and `beta0_white` the number of vacant matching-lattice components.
Let `q` be the common wrapping-difference event already used by P34
(`either`, and configuration-identically `cross`/`both`/`direction_*`).

Euler–Poincaré on this complex is the rewrite

```text
chi = V - E + F0 = beta0_black - beta0_white - q.
```

This is algebraically the same statement as the committed P34 identity

```text
C_black - C_white = q + V - E + F0.
```

Averaging recovers the expected Mertens–Ziff relation. The torus Euler
characteristic `chi(T^2)=0` is compatible with both empty and full
configurations, which have `chi=0`.

## What `q` is, and is not

`q` is a wrapping-**event** difference in `{-1,0,+1}`. It is not the homology
rank difference `r_black - r_white`.

The empty configuration is the counterexample on every named tiny quotient:

```text
q = -1,  r_black = 0,  r_white = 2,  beta0_black = 0,  beta0_white = 1.
```

The Euler identity still holds: `0 = 0 - 1 - (-1)`. Locked exhaustive counts
of the minority of configurations where `q` happens to equal `r_black-r_white`:

```text
axis L=2          4 / 16
axis L=3        162 / 512
gaussian-2-1     10 / 32
diamond L=2      68 / 256
```

## Cyclomatic bound, not a new observable

The graph-theoretic cyclomatic numbers

```text
kappa_black = E_primal - V_black + beta0_black
kappa_white = E_matching_white - V_white + beta0_white
```

satisfy `kappa >= wrapping rank` on every enumerated configuration. Wrapping
rank is a lower bound on the cycle-space dimension; it does not supply an
independent additive control.

## Variance-reduction branch

Issue #111 asked for exact control variates derived from Betti statistics, and
said to close that branch if they algebraically reduce to existing matching /
wrapping / motif controls. They do: every quantity above is a function of
`(V, E, F0, q, C_black, C_white)` already tracked by P34. No production
Newman–Ziff Betti accumulator is added.

## What this does not establish

- a new independent control variate or covariance identity beyond P34;
- a continuum homology interpretation or LCFT identification of `M`;
- any statement about Gaussian orientation amplitudes.
