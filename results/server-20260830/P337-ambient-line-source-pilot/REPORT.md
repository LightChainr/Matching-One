# P337 ambient-H1 source-column pilot

The new current-state coordinate is real, but the frozen observer pair does
not see it as a second source direction at 20k per size.

At fixed pre-insertion occupancy, the new column is

\[
W_{\rm line}(k)=\mathbf 1\{K_1\le k<K_2\}
\exp\!\left(4i\arg(P\ell)\right).
\]

It records the lifted primitive ambient-H1 line throughout the rank-one
plateau and is exactly zero off that plateau.  This is neither `qJ`, `q^2J`
nor a local motif.  The same permutation batches contain its products with
`O_far`, `O_sep4` and `JS`, plus `K1/K2/mark12_h4` for the independent
completion-cause readout.

## Frozen primary result

| size | determinant Re | determinant Im | chi2(2) | upper-tail p | normalized wedge |
|---|---:|---:|---:|---:|---:|
| N325 | -0.004528 +/- 0.006448 | 0.000588 +/- 0.009490 | 0.619 | 0.734 | 0.0231 +/- 0.0671 |
| N425 | 0.015843 +/- 0.015852 | -0.009578 +/- 0.011349 | 1.666 | 0.435 | 0.01638 +/- 0.01397 |

Neither size crosses the frozen alpha=0.01 `chi2(2)>9.21034` gate.  The two
sizes are independent evidence blocks and were not pooled.

This null is not caused by `W_line` being a disguised `JS`: after centered
same-batch Gram projection onto `JS`, `W_line` retains 98.27%/98.47% of its
norm at the two N325 orientations and 98.71%/98.95% at N425.  The pilot has
therefore separated two questions cleanly: the new state coordinate is almost
orthogonal to `JS` in configuration space, while its couplings to the present
`O_far/O_sep4` rows do not lift their source rank.

## Completion-cause control

The `mark12_h4` completion hazard was standardized on common
`(ell_u,ell_v,plateau age)` support.  Common-support coverage exceeds 99.5%
on every orientation.

| size | conditional mark12-H4 | SE | t(19) | risk-composition remainder |
|---|---:|---:|---:|---:|
| N325 | -0.000132 | 0.002798 | -0.047 | 0.000416 +/- 0.000653 |
| N425 | -0.000891 | 0.001917 | -0.465 | 0.000014 +/- 0.000300 |

Neither the conditional marked hazard nor the risk-composition remainder is
resolved.  Thus this pilot does not connect the finite-size completion H4 to
the new state column.

## Mechanism decision

Do not upscale the identical `O_far/O_sep4` lane from this pilot.  The useful
next move is a genuinely line-typed observer row (seam/winding incidence or a
landing observable addressed relative to the current lifted line), because
the missing rank is now localized to observer coupling rather than absence of
a new state coordinate.  This is a finite pilot statement, not a field or
exponent identification.

Both runs passed all endpoint, site, line, local-mark, index and separated-mark
exact audits with zero failures.  N325 ran on `TgFr7R` in 6.25 seconds and
N425 on `XPk2PZ` in 8.97 seconds.
