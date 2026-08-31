# P334/P429: cooperative branching closure from the same checkpoint

This is a fresh dynamic continuation experiment, not another static geometry
proxy. Its parent is P429 production commit `751f8b3`.

For a rank-one checkpoint C with d vacant sites, let b1 be the number of safe
one-site insertions and b2 the number of safe unordered two-site sets. For each
safe first site v let c_v count safe second sites. Monotonicity gives the exact
incidence identity `sum_v c_v=2*b2`. Set `s1=b1/d`,
`s2=2*b2/[d(d-1)]`, and `q_v=c_v/(d-1)`.

One common uniform update followed by two independent one-site clones has
conditional success `sum_v c_v^2/[d(d-1)^2]`. Uniform-safe-successor closure
would instead give `s2^2/s1`. Their exact difference is

```
Delta(C) = [b1 * sum_v c_v^2 - (2*b2)^2] / [d*b1*(d-1)^2]
         = s1 * Var(q_v | safe v, C) >= 0.
```

Set Delta(C)=0 when b1=0. The numerator is evaluated in exact integer
arithmetic, so positivity is not floating-point subtraction. The b2 census
already requires all c_v; saving their squared sum has no extra topology cost.
The primary ensemble estimand is the average of Delta(C), conditional on
rank-one C. This Rao–Blackwellizes the sampled common-update statistic without
changing its expectation.

Every row also retains H2, b2, the selected exact q_after and q_after^2, and
Y1/Y2/Y1Y2 from the two independently tagged clones. Their calibration rows
are evaluated jointly. A base permutation (including both paired orientations)
is one checkpoint cluster; the clones are never two sample units.

Pre-data decisions are in `analysis/p334_cooperative_closure_freeze.json`:
N325/Tg and N425/XP, 20k base permutations each, nominal alpha .01 two-coordinate
Wald gate with both coordinates positive, and no automatic extension.

The narrow eliminated mechanism, if positive, is uniform safe-successor
continuation at these finite checkpoints. It is not a proof of path memory,
nonlinear H2/b2 insufficiency, a scale law, or a continuum field identification.

## Result

Pending the frozen fresh production. No earlier P429 clone-gap sample is pooled.
