## #437 / #419: a genuinely readable fixed-support high-order topology endpoint

Branch `experiment/p437-fixed-support-lower-bound-20260831`; freeze/execution
`3c3fc57`. New raw block, full batch value classes/covariance and report:
`results/local-20260831/P437-N112-fixed-S5-20k/`.

**Main result:** for the previously fixed five bonds
`S={0,28,56,84,112}`, fresh 20k square-bond p=.5 backgrounds give

`B_S=E|D_SF|² = (3.23893 +/- .21639)e-6`, mean/SE **14.97**.

There are **384 nonzero backgrounds, 1.92%, all 100 batches represented**.
All 32 support states are enumerated in each background. No support scan,
importance-combination factor, or old-stream extension is involved.

The normalized derivative is `D_S=2^-5 Delta_S`. Exactly,
`B_S=sum_{T contains S}|Fhat(T)|²`, and the original population high-pass
energy obeys **`A_HP >= (9765/32768) B_S`**. The measured right-side parameter
is `(9.65215 +/- .64485)e-7`; this estimate is **not** a statistically certain
numerical lower bound.

**Mechanism space changed:** a collective fifth-order local support has
non-negligible, directly readable topology response in the true product
distribution, beyond the exceptional exact witness. It includes all spectral
orders containing S, not a pure degree-five field or total high-pass energy.

**Typed localization:** child0 supplies 86.1% of sampled energy; the support
itself is not child-symmetric, so this is a fixed-support local response—not
a global C3 preference. The 20k block contains no simultaneous multi-child
activation; this is observed, not an exact selection rule.

**Cost/dependency/stop:** 36.82 local wall seconds, 333.14 CPU seconds; new group
`p437-N112-fixed-S5-lower-bound-fresh20k-20260831`. Old six-noise stop remains
unchanged. No new Huawei job, support rotation, expanded sample, PR or merge.
