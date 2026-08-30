# Finite Bernoulli model-confidence e-process

This Issue #126 slice inverts two exact likelihood-ratio e-processes into a
confidence sequence over the predeclared model family `{p=1/3,p=2/3}`.  Each
model is retained at time `t` when the other-versus-null likelihood ratio is
below `1/alpha=4`.

The four-step Fraction oracle enumerates every binary path under each possible
declared truth.  It verifies that each e-value has expectation one at times
0 through 4, that the e-value stopped at its first threshold crossing (or the
horizon) still has expectation one, and that confidence-set exclusion occurs
exactly on the crossing paths.  Under either true model,

```text
P(the true model is excluded at least once) = 13/81 < 1/4.
```

Thus the simultaneous coverage over the four declared times is at least 3/4.
A model may re-enter because the confidence set uses the current e-value; the
guarantee concerns never excluding the truth over the whole sequence.

Run:

```text
python3 scripts/bernoulli_model_confidence_eprocess.py
python3 -m unittest tests/test_bernoulli_model_confidence_eprocess.py -v
```

## Boundary

Coverage is conditional on one of the two predeclared models being true.  There
is no coverage statement when data come from an off-model distribution, and
the oracle does not permit generating new candidate models adaptively from the
same observations.  It chooses no production allocation or scientific model
and does not complete Issue #126.
