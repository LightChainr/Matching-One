# Hexagonal CM degree-2 Hecke fingerprint

Issue #164 asks for geometry fingerprints that separate a typed thermal
spin channel from generic scalar or Jordan contamination.  At the hexagonal
CM point this can be frozen before any Monte Carlo run.

We use

\[
 \omega=e^{i\pi/3},\qquad \zeta=e^{2\pi i/3}=\omega^2,
 \qquad \widehat E_k(\tau)=\operatorname{Im}(\tau)^{k/2}E_k(\tau),
\]

with period basis `(1,tau)` and child order

\[
 2\omega,\quad \omega/2,\quad (\omega+1)/2.
\]

Direct high-precision q-series evaluation, the modular transformation law,
and the normalized degree-2 Hecke equation independently give

\[
 (\widehat E_4(2\omega),\widehat E_4(\omega/2),
   \widehat E_4((\omega+1)/2))=A(1,\zeta,\zeta^2),
\]

while `Ehat_4(omega)=0`, and

\[
 \widehat E_6(2\omega)=\widehat E_6(\omega/2)
 =\widehat E_6((\omega+1)/2)=\frac{11}{4}\widehat E_6(\omega).
\]

For E6 the corresponding holomorphic ratios are `11/32, 22, 22`.
The q-series script adaptively truncates each sum and records the terms used;
the frozen JSON retains 70 significant digits.

The E4 triplet has reflection-even real pattern `1:-1/2:-1/2` and
reflection-odd imaginary pattern `0:+sqrt(3)/2:-sqrt(3)/2` in this basis.
A child difference or the nontrivial three-point DFT therefore cancels any
common scalar child contribution.  The E6 triplet is instead a constant
vector.  Its `11/4` child/parent score is meaningful only after the scalar
normalization requested in Issue #161; equal children alone cannot distinguish
spin 6 from scalar contamination.

Finally, `11/4` is a **shape ratio** after area normalization.  It must not be
substituted for similarity scaling: under a scalar similarity by `m`, a
weight-`k` quantity still scales as `m^{-k}`.  These are exact modular-form
identities, not by themselves an identification of a percolation observable.

Reproduce with:

```bash
python scripts/derive_hexagonal_degree2_hecke.py --dps 90
python -m unittest tests.test_hexagonal_degree2_hecke
```
