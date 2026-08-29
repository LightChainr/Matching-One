# Issue #200 Phase A: N580 production runner

The N580 designs `(24,2)` and `(18,16)` have gcd two.  Their period matrices

```text
((24,-2),(2,24)), ((18,-16),(16,18))
```

both have Smith invariants `(2,290)`.  They therefore cannot be added to the
cyclic Gaussian quotient implementation in `threshold_rank_orientation_mc.cpp`.
The schema-compatible general integer-period backend is the exact existing
implementation for these graphs.  This branch adds N580 as a built-in design
there; the production command names the compiled executable
`threshold_rank_orientation_mc` so the established histogram/moment analysis
chain consumes it without conversion.

The exact C++ self-test now locks both N580 matrices and their Smith class.
The Python regression independently constructs each period quotient, compares
the C++ threshold histograms for both orientations, and verifies byte-identical
integer outputs with one and two workers.

The 100M/100-batch/8-thread contract is
`experiments/p200_n580_phaseA_100m_20260829.yaml`.  Seed `2026102001` and counter
interval `[12000000000,12100000000)` are outside the P50 and P154 stream
domains.  The two N580 orientations share each counter-derived HNF-label
permutation; worker scheduling cannot change output rows.

The downstream score input is frozen separately in
`predictions/p200_n580_q2_jordan_score_input_20260829.json`.  It preserves the
four-state order `(I_S,I_Du,T_D,T_Su)`, the full propagated q2/Jordan source
covariance, and the rule that N580 is one joint four-coordinate score.  The raw
runner does not inspect targets and cannot select a model.

This Phase A branch intentionally contains no N650 path flags.
