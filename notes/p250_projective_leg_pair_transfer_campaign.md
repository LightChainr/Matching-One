# P250 projective-leg pair-transfer campaign

## Decision

The fresh 10k result is more useful as a discovery of a **mesoscopic charged
propagator** than as a failed cubic/OPE experiment.  This campaign therefore
drops the cubic row completely and records the complex two-point transfer row
at separations 1--6 with one joint batch covariance.

For hand `h` and charge `r=1,2`, define

\[
T_{h,r}(d)=\frac12\left\langle
O_{h,r}(0)O_{h,-r}(d\hat x)+O_{h,r}(0)O_{h,-r}(d\hat y)
\right\rangle .
\]

The existing Hermitian denominator `G_r(d)` is exactly `Re T_r(d)`, so the
fresh 10k stream can select a sample size without looking at a new transfer
shape.  Its weakest real-row resolutions at `d=1,2,3` are respectively
15.397, 8.608, and 2.504 standard errors.  The first frozen grid point whose
square-root forecast clears five standard errors at all three distances is
40k.  Separations 4--6 are informative tail diagnostics, not sample-size
requirements.

## Frozen amplitude-free decision

Let `R12=T(2)/T(1)` and `R23=T(3)/T(2)`.  The primary score uses all four
hand-charge channels and their joint delete-one-batch covariance.

- A single exponential transfer eigenvalue predicts
  `m12=-log|R12| = m23=-log|R23|` and a constant complex phase step.
- A power tail predicts
  `eta12=m12/log(2) = eta23=m23/log(3/2)` and a constant phase step.
- The complex arguments separately test whether a deck-character phase is
  nonzero and whether its adjacent step is constant.

This separation is deliberately sharper than a fit with a free amplitude:
exponential scale, power-like tail, and deck phase can fail independently.
If both simple shape scores fail while the pair row remains resolved, the
result points to a finite-torus transfer-state mixture rather than to a null
observable.

## Scientific card

- **Mechanism changed:** the projective-leg insertion becomes a propagating
  charged state rather than an OPE-support proxy.
- **Not proved:** locality, primary-field identity, or a universal OPE
  coefficient.
- **Observer / sector / source / geometry:** complex `Z5` pair row; charges
  1 and 2; common Bernoulli counter stream; both projective-leg hands on the
  N325 parent torus.
- **Dependency group:** the fresh stream is independent; its 40k design uses
  only denominator variance from the frozen 10k archive.
- **Next upweighting signal:** one of the amplitude-free shape residuals closes
  jointly, or a stable nonzero deck phase survives both charges and hands.
