# Aggregate TM is a curvature-corrected Rayleigh inequality

Choose a uniform fixed-line state and an ordered pair of distinct absent sites. The two-site face has four possible motifs: coexit `D`, mixed `M`, synergy `Y`, or flat `F`; write `T=D+M+Y+F`.

Direct integer expansion gives

`M^2 + 4Y(T-D) - 4DF = 4(m-1) TM_margin`.

Equivalently, with `p` the one-mark exit probability, TM is `P(D) <= p^2+P(Y)`. The ordinary Rayleigh determinant is corrected by exactly the concave two-site faces where neither single exits but the double insertion does.

## The two covering mechanisms are both essential

Ordinary Rayleigh `4DF<=M^2` fails 16 rows. Synergy-only coverage fails 68 rows. Their corrected sum passes all 984 rows and satisfies the integer identity in all 984 rows. The finite motif cone has 9 Pareto-minimal cover rays; no single cover term spans it.

## Why aggregation over every displacement is essential

The corrected inequality fails on 3900 of 51912 individual site-pair tables, beginning at N=6. Even grouping by quotient order leaves 220 failures in 3330 order tables. Thus delta-by-delta and inversion/order pairing are too strong.

A Fourier sum-of-squares route is also false: an exact two-site contrast is negative in 802 rows. Alexander reflection does not identify the two margins: only 18 of 492 reflected primal/matching pairs are equal, because an exit face reflects to a birth face.

## Single remaining motif inequality

The general topology target is now one explicit four-face injection:

`coexit x flat -> mixed x mixed  OR  synergy x non-coexit`.

Proving its aggregate cardinality `4DF<=M^2+4Y(T-D)` after summing all relative displacements is exactly aggregate TM. By the regular-cut theorem, it then supplies the one-mark Hall injection and every proper Hall cut automatically.
