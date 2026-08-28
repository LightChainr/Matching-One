# Tiny exact matching-polynomial complex-zero pilot

The map uses exact enumeration coefficients and 100/130-digit numerical roots.
It is an algebra-and-design pilot, not a thermodynamic zero-density claim.

## Root audit

| geometry | L | N | degree | real/nonreal | physical root | imag RMS | max residual | max conjugate error | max matching error |
|:---|---:|---:|---:|:---:|---:|---:|---:|---:|---:|
| axis | 1 | 1 | 1 | 1/0 | 0.5 | 0.0 | 0.0 | 0.0 | 0.0 |
| axis | 2 | 4 | 4 | 4/0 | 0.541196100146 | 0.0 | 8.8e-121 | 0.0 | 2.0e-120 |
| axis | 3 | 9 | 9 | 3/6 | 0.586511455113 | 0.40871515 | 1.29e-120 | 0.0 | 4.8e-120 |
| axis | 4 | 16 | 16 | 4/12 | 0.590672112331 | 0.4938347 | 3.56e-120 | 0.0 | 5.0e-120 |
| diamond | 1 | 2 | 2 | 2/0 | 0.707106781187 | 0.0 | 6.54e-121 | 0.0 | 2.0e-120 |
| diamond | 2 | 8 | 8 | 4/4 | 0.604563277854 | 0.30435964 | 2.87e-120 | 0.0 | 4.0e-120 |
| diamond | 3 | 18 | 18 | 6/12 | 0.594252321169 | 0.37479069 | 1.79e-120 | 0.0 | 4.0e-120 |

The matching column compares roots of `P(p)` against independently solved roots of
`-P(1-p)`. A single axis/diamond polynomial is not assumed to be self-matching.

## Declared train/holdout diagnostics

Each target is predicted from the preceding two sizes with `a+b/N`. This deliberately
cheap rule is scored before any richer fit; the exercise defines quantities that a future
zero-map computation can falsify instead of merely producing a visually suggestive cloud.

| geometry | metric | train L | held-out L | prediction | observed | signed error |
|:---|:---|:---:|---:|---:|---:|---:|
| axis | physical_root_0_1 | 2,3 | 4 | 0.602371829351 | 0.590672112331 | -0.011699717 |
| axis | imaginary_rms | 2,3 | 4 | 0.551765457373 | 0.493834697369 | -0.05793076 |
| axis | nonreal_fraction | 2,3 | 4 | 0.9 | 0.75 | -0.15 |
| diamond | physical_root_0_1 | 1,2 | 3 | 0.585573740199 | 0.594252321169 | 0.008678581 |
| diamond | imaginary_rms | 1,2 | 3 | 0.360722538168 | 0.374790685736 | 0.014068148 |
| diamond | nonreal_fraction | 1,2 | 3 | 0.592592592593 | 0.666666666667 | 0.074074074 |

## Frozen next-size predictions

| geometry | target L/N | metric | prediction |
|:---|:---:|:---|---:|
| axis | 5/25 | physical_root_0_1 | 0.59259790224352 |
| axis | 5/25 | imaginary_rms | 0.53323288619515 |
| axis | 5/25 | nonreal_fraction | 0.78857142857143 |
| diamond | 4/32 | physical_root_0_1 | 0.59064348632884 |
| diamond | 4/32 | imaginary_rms | 0.39944155119135 |
| diamond | 4/32 | nonreal_fraction | 0.725 |

Axis L=5 requires 2^25 configurations and remains within the reference engine's hard
limit; diamond L=4 has N=32 and needs a frontier exact engine. Preserve these predictions
unchanged if either target is later computed.

## What this pilot says

- Every solved polynomial passes conjugate pairing, dual matching-root pairing and
  independent-precision stability by many orders of magnitude.
- The unique real root in `(0,1)` rapidly approaches the known threshold from opposite
  sides for the two orientations; the two-point tiny-size extrapolations are visibly biased.
- Complex roots already proliferate by axis L=3 and diamond L=2, but the cloud summaries
  are not yet smooth enough to justify a conformal or Lee-Yang interpretation.
- The next useful step is one new exact size, scored against the frozen scalar targets above,
  before adding plots, clustering rules or modular interpretations.
