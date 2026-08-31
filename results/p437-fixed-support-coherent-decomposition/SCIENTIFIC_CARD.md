## #437 / #419: the 15-SE local energy is mostly degree>=6, with a stronger same-data bound

Zero new samples: reaggregate the signed classes of `386db0a`, fixed
`S={0,28,56,84,112}`. Branch
`analysis/p437-fixed-support-coherent-decomposition-20260831`.

- **Spectrum:** `mu=E[D_SF]=Fhat(S)` is the single exact-degree5 coefficient;
  `B_S-|mu|²` contains only strictly larger supports, hence degree>=6.
- **Measured coefficient:** `mu=(-7.10938e-5)+i(1.35316e-6)`, SEs
  `(1.13474e-5,4.36961e-6)`, full 2x2 covariance saved. Real component is
  -6.27 SE; the imaginary component is not used to infer a symmetry.
- **Exact phase:** reflection x->-x plus translation (9,0) maps child1 and
  its fixed support to child2. Thus their real coefficients are equal and
  **Im mu=0 exactly**. Real mu is permitted; the exact transport enumeration
  finds no anti-invariant support map forcing the whole coefficient to zero.
- **Energy decomposition:** the cross-batch U-statistic subtracts the
  squared-mean bias. Coherent weight is `(4.9083 +/- 1.6162)e-9`;
  outside-dependent weight is `(3.23402 +/- .21641)e-6`.
  The latter accounts for **99.8485% +/- .0511 percentage points** of B_S.
- **Stronger population inequality:**
  `A_HP>=h5|mu|²+h6(B_S-|mu|²)`, with exact `h6/h5=63/32`.
  Its RHS parameter estimate is **`(1.89885 +/- .12696)e-6`**, about
  **1.96728 times** the preceding h5 B_S estimate. Uncertainty uses the full
  joint covariance. This is not a statistically certain numerical lower bound.

This is secondary reuse of the same dependency group, not independent evidence
or a change to the frozen B_S primary. It explains why a local five-bond
insertion can have a strong topology response while most of its energy lives
in larger collective supports. No new MC, rotated support, field name or PR.
