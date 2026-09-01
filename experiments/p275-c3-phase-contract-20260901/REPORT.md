# P275 C3 phase-identifiability gate

Status: `CONDITIONAL_DESIGN_SOLVED__OBSERVABLE_TRANSPORT_MISSING`

This is an exact forward-design calculation for the same phase-calibrated
real C3 observable at two physical rotations.  It does not fit an empirical
field identity or authorize new sampling.

For a rotation difference `delta`, the real H4 and aliased H8 design spaces
have joint rank four exactly when

```text
sin(6 delta) != 0.
```

Under a shared unknown complex amplitude, the fraction of an isotropically
balanced wrong-model signal left after profiling the other model is
`sin(6 delta)^2`.  This immediately separates three amplitude contracts:

| Transport contract | Exact conclusion |
|---|---|
| same complex amplitude | any `delta` not divisible by 30 degrees is identifiable |
| unknown nonzero signed real gain | 7.5 degrees separates H4/H8; 15 degrees can alias after a sign reversal |
| arbitrary complex gain | no two-rotation design can identify H4 versus H8 |

The 7.5-degree design has wrong-model residual fraction `1/2` under the shared
amplitude contract and signed-gain phase separation `sin(12 delta)^2=1`.
The 15-degree design has residual fraction one for a fixed shared amplitude,
but the exact counterexample

```text
amplitude 1+2i, H4 gain +1, H8 gain -1
```

produces the same two readouts.  Full-rank covariance whitening changes the
amount of separation but not these rank statements.

The current P0 is therefore narrower than “measure another angle.”  Theory
must first map the original observable, its normalizer and pooled moving root
to one of the three amplitude contracts.  If arbitrary complex gain remains
allowed, this C3 route is formally unidentifiable and should be downgraded.
If shared amplitude or signed real gain is proved, freeze the corresponding
two-vector prediction and score existing source-matched assets before any
new acquisition.

`RESULT.json` records 66 exact symbolic checks, including 12 direct C3
Fourier checks, rational Gaussian-multiplier phase arithmetic, alias
counterexamples and covariance-weighted controls.  It was generated with:

```text
/Users/lc/python-envs/research-py311/bin/python phase_design.py --out RESULT.json
```

No lattice enumeration, Monte Carlo block, geometry choice or empirical C3
reanalysis was performed.
