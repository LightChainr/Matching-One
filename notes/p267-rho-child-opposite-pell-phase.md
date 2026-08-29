# Opposite-Pell N60 completion-phase discriminator

The N112 reveal establishes r0+r1 support but leaves a coherent r2-shaped
failure of the particular E4+E6 functional relation.  The next experiment
therefore tests that missing complex relation, not the three pure rows again.

Use all three degree-two children of the opposite-side Pell parent
`[[6,3],[0,5]]`:

```text
2*tau       [[6,6],[0,10]]
tau/2       [[12,3],[0,5]]
(tau+1)/2   [[12,9],[0,5]]
```

For each three-child vector define the exact complex covector `w` by

```text
w dot E4 = w dot E6 = 0,
c_N = (w dot delta_H4)/(w dot E4^2).
```

This is normalization free up to a positive real radial amplitude.  The
primary E4-squared completion is the positive ray `c60=a*c112`, `a>0`.
The signed-Pell opponent is the negative ray `c60=-a*c112`.  Each ray profiles
one positive scale and a latent complex reference coefficient against the
joint N112/N60 covariance, leaving one residual degree of freedom.  Both
rejected means phase incoherence; both surviving means insufficient phase
resolution.  The decision threshold is frozen at `.01`.

The N60 and N112 annihilators both contract E4-squared with phase `pi/2`
and magnitudes `729.49046` and `731.72966`, respectively.  Thus phase
preservation is an actual prospective prediction, not a post-reveal choice.
No pilot coefficient, amplitude, exponent or phase is fitted before the
two-million-sample acquisition.

The full primitive-line sum passes the exact scalar/spin-4 alias gate at all
three children.  Tiny exhaustive classification has zero invariant failures.
The production manifest was authorized only after implementation commit
`0a25a46ab34a662e6561c4b2cc9a3d0c21bffbd1` recorded the runner, scorer,
tests and frozen CLI/counter fields.  This authorization changes no model or
acquisition parameter.

## Frozen reveal

The six-million-sample Huawei acquisition gives

```text
c112 = -5.725761e-5 + 2.394937e-6 i
c60  = -3.602105e-5 + 9.831777e-6 i
c60/c112 = 0.635176 - 0.145144 i
```

The raw phase difference is `-0.22465` radians.  The phase-preserving ray
has an interior positive scale `.79676` and survives (`chi2=.95841/1`,
`p=.32759`).  The signed-Pell flipping ray has `p=.03323`, above the frozen
`.01` rejection threshold, but its best scale is the lower numerical boundary
`1e-6`; it survives only by collapsing to the unresolved N60-zero model
(`p=.10363`).

Therefore the preregistered decision is `phase_unresolved_both_survive`.
E4-squared phase preservation remains viable and is the only model with an
interior nonzero amplitude.  A sign-flipped or Jordan/incoherent completion
is not selected, but cannot be excluded at the frozen threshold because the
N60 annihilator coefficient itself is not resolved.  This is an informative
power boundary, not permission to rescore at `.05`.
