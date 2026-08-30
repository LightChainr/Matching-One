# Preregistration: N340 same-lineage second child

Frozen before generating or inspecting N340 data. Multiplication of N170 `(11+7i,13+i)` by `1+i`, followed only by D4 canonicalization, gives N340 `(18+4i,14+12i)`. The exact H4 covector flips again; the projective scalar remains a frozen zero-amplitude control.

Primary coordinate: `A_H=(K_second-K_first)/(c_second-c_first)`. Fixed targets:

- `nominal_area_H4`: `A_H=-0.0036024562797`, pair `+0.00574398565289`
- `observed_N85_to_N170_effective`: `A_H=-0.00769210003653`, pair `+0.0122647740375`
- `scale_neutral`: `A_H=-0.0111114941452`, pair `+0.0177168737098`

The secondary effective transfer is the already revealed N85-to-N170 ratio `0.692265139`; it is not an exponent fit. N85 target uncertainty is retained.

Design: 12M/shape, 80 aligned batches, seed `202608337340`, HZsCM6. N170 variance projects H4-amplitude SE `0.0007966`; fixed target gaps are 5.13, 4.29, and 9.43 measurement SE.

Scoring order is frozen: nominal area H4, observed effective continuation, scale-neutral. Both measurement-only and source-uncertainty-aware residuals are reported. No exponent fit, basis change, or H4/H8 revote is allowed.
