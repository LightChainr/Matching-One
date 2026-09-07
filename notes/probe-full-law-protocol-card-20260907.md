# Protocol card — full-law affine gauge (#617, P7 freeze), 2026-09-07

```text
declared chart   = span{1, Q_base, g_frozen}
g_frozen         = #582/#584 consensus, PR #614 commit f0981a98
                   (results/p582-amplitude-law/latest.json; never re-fit)
transport        = #612 identity, applied BEFORE any residual:
                   a_curv = 2/(h0+h1) [a1 - a0/k + ell_C(r1) - ell_C(r0)/k],
                   k = 1 + h0*beta0; general attachment via the covector
                   expansion (exact, 1.6e-13)
weightings       = spin0 AND equal, always together (paired reporting)
reconstruction   = production inverse-CDF: nine deciles u in {0.1..0.9},
                   binomial-profile CDF (P5: this readout is load-bearing;
                   the 4.25% lives only here)
residual rule    = sign-flip label null (10^4 draws, seed 6170425) on any
                   claimed shape; P7 verdict on the 4.25%: DIES (p=0.647)
forbidden        = second exponent; chart change after seeing a residual;
                   folding 725 into a three-size curvature
surviving
invariant        = I4 only: projective rho = delta/A4, predeclared
                   denominator, attachment-stable (15.3% span < 50%, P10),
                   weighting-stable (0.02%, P8)
forecast standing = N=725 finite forecast SUPPORTED in both weightings
                   (3-se interval inside the pre-registered +/-5% band);
                   this does not make the model exact
holes            = N-normalisation axis unswept; weighting fibre 1-D
                   (no orthogonal direction exists); p50 two-size (no
                   curvature); direction transport across readout ranks
                   undetermined (grid-bound g)
status sentence  = see notes/probe-full-law-affine-gauge-20260907.md (P11)
```
