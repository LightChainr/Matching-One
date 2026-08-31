# Fixed Q1 trace scorer: code staged before theory-gate release

This entry implements only contract `964ef2032effbe59f9158c158cf06a2c0844d7ee`.
**Do not execute it until the packing and generic-character proofs are
both complete and the coordinator explicitly releases the calculation.**

The two fixed coefficient packets on the already complete mod6 histograms
are `beta1=-1_A-1_B` and
`H=(K+g+3)/2*1_A+(K+g+1)/2*1_B`, with A/B exactly as the contract names
them. The code reports any unsupported `(q,bad2,n_bad3)` row and stops
before scoring, without changing the source.

It reuses the existing exact interval geometry evaluator and three-term
normalization map. The baseline coefficients now have **iid** weights
`count*y^K`, not the old Q4 `count*2^-g*y^K`. The already published Q1
p-root is converted once to `y=p/(1-p)` and the published thermal slope
to `D_y=D_p/(1+y)^2`; the published U/A is unchanged. No root finder or
Sstar/Bvac response routine is invoked.

The primary is an epsilon-insertion response at Q1. The secondary is an
additive raw-trace derivative attribution in the fixed reduced partition
`y^K Q^(-(K+g)/2)`, not an invariant share and not the mixed derivative
`d_Q d_epsilon U`. Their common data do not supply independent evidence.

After release, substitute the two supplied proof commits:

```bash
/Users/lc/python-envs/research-py311/bin/python \
  scripts/p337_q1_trace_continuation_score.py \
  --packing-gate-commit PACKING_SHA \
  --character-gate-commit CHARACTER_SHA \
  --output-dir results/p337-q1-trace-continuation
```

Before gate release, only Python AST syntax parsing is authorized for this
new script. No counts or new response are consumed by that syntax check.
