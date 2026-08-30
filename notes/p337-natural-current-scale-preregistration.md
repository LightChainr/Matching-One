# Preregistration: natural charged-current scale test

This artifact is frozen before generating or inspecting N85 projective-birth data.

```text
K_A = p(1-p) Jminus_A/W_A = d_eta log W_A,
Delta_K_A = K_A(second)-K_A(first).
```

N65 gives `K_first=-0.0415636809311`, `K_second=0.0279169787293`, and `Delta_K=0.0694806596604 +/- 0.0204` by aligned batch jackknife.

Frozen N85 targets:

- zero: `0`
- source-fitted scale-neutral: `0.0694806596604`
- source-fitted project H4, `(85/65)^(-13/8)`: `0.0449306170854`

No old N85 file contains the required `tau1,ell,tau2` statistics. The fixed next pair is `(9+2i,7+6i)`, 200000 samples per shape, 20 aligned batches, seed `202608337`.

N65 variance projects an N85 SE of `0.006441`. The 200k block exceeds the 3-sigma scale-neutral versus H4 gap requirement while remaining a short fresh run.

Scoring will report both residuals to the frozen numeric targets and predictive residuals including N65 fitted-target variance. N85 cannot choose a new normalization, exponent, or target.
