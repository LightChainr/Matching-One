# Gaussian doubling test: a parameter-free spin-4/exponent check

A useful held-out test follows from multiplying both Gaussian period generators by `1+i`.

For a primitive orientation represented by `z=a+ib`, multiplication by `1+i` gives

\[
z'=(a-b)+i(a+b),
\]

so the norm doubles:

\[
N'=|z'|^2=2|z|^2=2N,
\]

and the physical orientation rotates by `pi/4`.

For a pure spin-4 correction

\[
\Delta M_N=C\,\Delta\cos(4\theta)\,N^{-13/8},
\]

this operation flips the angular harmonic exactly:

\[
\Delta\cos[4(\theta+\pi/4)]=-\Delta\cos(4\theta).
\]

Therefore the `x=21/4` hypothesis predicts the parameter-free ratio

\[
\boxed{
\frac{\Delta M_{2N}}{\Delta M_N}
=-2^{-13/8}
\approx-0.3242098886627524.
}
\]

No fitted amplitude and no value of `p_c` enter this ratio.

## Two immediate lineage pairs

The current same-N designs have exact doubled descendants.

### N=65 -> N=130

Original representations:

- `(8,1)`
- `(7,4)`

After multiplying by `1+i`, canonical representatives are

- `(9,7)` corresponding to `(8,1)`;
- `(11,3)` corresponding to `(7,4)`.

Keep this lineage ordering when forming `Delta M`, even if the display order is reversed.

Using the current `N=65` production value `Delta M ~= +1.00377e-3`, the no-fit prediction is approximately

\[
\Delta M_{130}\approx-3.25\times10^{-4}.
\]

### N=85 -> N=170

Original:

- `(9,2)`
- `(7,6)`.

Doubled/rotated descendants are

- `(11,7)` corresponding to `(9,2)`;
- `(13,1)` corresponding to `(7,6)`.

The current `N=85` value `Delta M ~= +7.6033e-4` predicts

\[
\Delta M_{170}\approx-2.47\times10^{-4}.
\]

## Why this test is unusually strong

A generic power-law fit can trade exponent against amplitude and subleading corrections. This construction fixes both the scale ratio and the spin phase:

- scale factor in physical length: `sqrt(2)`;
- spin-4 phase: exactly `-1`;
- predicted magnitude ratio: exactly `2^-13/8` under the leading hypothesis.

So one comparison tests the angular spin assignment and scaling dimension simultaneously.

## Protocol

1. Use the existing same-N fixed-p engine first; no new full-curve code is required.
2. Preserve exact lineage ordering so the predicted sign is not erased by canonical `(a>=b)` display conventions.
3. Use a fresh independent seed.
4. Do not fit an exponent to these points before evaluating the fixed prediction above.
5. Report the standardized residual of

\[
\Delta M_{2N}+2^{-13/8}\Delta M_N.
\]

6. After C05 exists, repeat the same test on the thermal-even spin-4 projected full-curve observable, which is the cleaner operator-level target.

## Interpretation

Agreement would be strong evidence for a spin-4 correction with `N^-13/8=L^-13/4` scaling. Failure would distinguish several possibilities: significant logarithmic/subleading corrections, the wrong exponent, more than one spin-4 sector, or an accidental amplitude collapse at N=65/85.
