# PR267 E_top tau x topology-map N50 factorial

The exact minimum four-cell design crosses `tau=i,2i` with
`Smith=(1,50),(5,10)` at fixed determinant 50. One rational rotation
`O=(1/5)[[4,-3],[3,4]]` maps cyclic to noncyclic in both rows.

| field | P4 at tau=i | P4 at tau=2i | interaction | interaction SE |
|---|---:|---:|---:|---:|
| A_top | 0.00130445 | 0.004064855 | 0.002760404 | 0.00093 |
| E_top | -0.0004205952 | -0.0009717124 | -0.0005511172 | 0.0007 |
| C | -0.000557946 | -0.00189451 | -0.001336564 | 0.00015 |
| W | 0.0002112082 | 0.001094431 | 0.0008832226 | 0.00026 |

## Frozen primary

Character-normalized interaction: `chi2=236.756/4`, `p=4.63371e-50`.
The 20k pilot crossed the frozen gate at `chi2=59.7457/4`; the displayed score uses the authorized
100k total per missing cell (20k pilot plus disjoint 80k extension).
The additive no-interaction factorial is therefore eliminated: the topology-map
response changes with tau at this fixed N.

The raw modulus and topology-map main-effect scores are retained in `score.json`,
but they are descriptive because their four coordinates are strongly correlated views
of the same threshold clocks. The primary endpoint is only the four-vector interaction.

## Fixed A+C diagnostic

Using the previously frozen coefficients gives residual `0.001617755 +/- 0.00078` (`z=2.075`, `p=0.038`). It survives
the existing alpha=.01 boundary but is mildly tense at .05. One interaction equation
cannot identify both A and C coefficients, so this is compatibility rather than pinning.

The P00 anchor histogram and moments reproduce byte-identically across Zy and XP
within both the pilot and extension blocks.
This is a finite-N geometry factorization result, not an exponent, root-character or
asymptotic modular-law claim.
