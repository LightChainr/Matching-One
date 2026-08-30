# P250 cross-scale projective-leg freeze

## Geometry choice

The N65 parent has side length `sqrt(65)=8.062`, so only `d=1,2,3` form a
clean local window; `d>=4` is already antipodal or periodic-image dominated.
An exact enumeration of primitive oblique Gaussian parents above norm 65
selects `10+i`, norm 101, as the first parent with five distinct signed
unit-axis separations before quotient identification.  Its norm-five children
have order 505 and periods `(19+12i)` and `(21-8i)`.

Norm 100 is intentionally not used: its real, nonprimitive parent makes the
two norm-five hands a reflection pair and changes the hand semantics.  N101
is therefore the smallest enlargement that preserves the current oblique
same-parent construction and gives `d=1..5` without displacement collisions.

## Models frozen before the N101 stream

All models are amplitude-free ratios to `T(1)`.  No N101 parameter is fitted.

1. **Thermal cylinder:**
   `[sin(pi d/L)/sin(pi/L)]^(-5/4)`.  This is the parameter-independent
   `2x=5/4`, `x=5/8` target requested by the thermal-family interpretation.
2. **Source-fitted cylinder:** the same sine law with one common exponent
   estimated solely from the old N65 `d=1..3` rows.  The frozen value is
   `2x=1.4145386844 +/- 0.0526456379`.
3. **Finite-correlation alternative:** one lattice-unit exponential whose mass
   is estimated solely from the same source rows,
   `m=0.7120974786 +/- 0.0267089600`.

The N101 target scores `d=2..5` jointly with the full hand-charge covariance;
`d=5` is also reported separately as the decisive antipodal holdout.  Source
fit uncertainty is propagated independently into the latter two predictions.
This separates a massless sine-law scale, a fixed lattice correlation length,
and failure caused by periodic images or a transfer-state mixture.

## Sample size

The old N65 40k batches supply only a variance proxy and the fixed-5/4 sample
design amplitude.  Using the largest per-replica real-row variance over
`d=1..3`, the predicted weakest N101 `d=1..5` resolution is 4.77 at 20k,
6.74 at 40k, and 9.54 at 80k.  The frozen choice is the first passing point,
40k.  The source-fitted exponent, exponential mass, and all phase information
were excluded from this choice.

## Scientific card

- **Mechanism changed:** tests whether the resolved charged propagator carries
  a geometry-scaled massless law or a fixed lattice correlation length.
- **Not proved:** a local thermal primary, a universal common exponent, or any
  cubic/OPE coefficient.
- **Observer / sector / source / geometry:** complex projective-leg pair row;
  Z5 charges 1/2 and both hands; independent fresh counter stream; N101 parent
  with N505 children.
- **Dependency:** N65 supplies frozen parameters and covariance only; N101 is
  an independent geometry and seed.
- **Upweighting signal:** the parameter-free 5/4 sine law predicts the joint
  N101 row and held-out d5 while the source-mass exponential fails.
