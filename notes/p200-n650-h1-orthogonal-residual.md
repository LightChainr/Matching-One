# N650 typed-H1 orthogonal residual

This is a post-reveal model-elimination score on the existing 20,000-sample
N650 mixed-join archive. It adds no Monte Carlo and does not rescore the N580
q2/Jordan clock. The runner has already subtracted the isolated-fibre
`C2 x C5` incidence background configurationwise:

\[
R_c=J_c^{\rm full}-J_c^{\rm local}.
\]

The previous context map established large common `ES/OS` residuals and null
`ED/OD`. The new question is narrower: can the surviving common signal be only
one endpoint ambient-H1 mechanism with a shared conversion gain?

## Gain-free H1 alignment test

Let `P=(primary_ES,primary_OS)` and
`H=(ambient_ES,ambient_OS)`. A one-gain endpoint-H1 explanation requires
`P=beta H`, hence the basis- and gain-free determinant

\[
P_O H_E-P_E H_O=0.
\]

The complete shared-batch 8x8 delete-one covariance of primary and ambient-H1
is retained. The observed determinant is

```text
68.95907262 +/- 0.49477963, z=139.37.
```

Equivalently, fixing `beta=primary_ES/ambient_ES=26.86669918` leaves

```text
H1-orthogonal OS = -34.49592187 +/- 0.24737557, z=-139.45.
```

Both even-to-odd split-half directions retain the same residual
(`-34.3618`, `-34.6304`). The one-gain ambient-H1 endpoint mechanism is therefore
eliminated. This is not merely the old large `OS`: it certifies that `OS` is
not the typed projection of the observed ambient-H1 common vector.

After applying the same frozen gain, the remaining geometry differences are

```text
ED residual = -0.00288833
OD residual = -0.00325164
joint chi2(2)=0.01976, p=0.99017.
```

Thus the new direction is common across the two static HNF embeddings and
matching-odd in the typed colour layer. A saturated model with unrelated even
and odd gains is not testable from only these two common rows and is not used
as an explanation.

## Scientific boundary

This changes the mechanism map: exact local incidence plus one ambient-H1
endpoint direction does not exhaust the N650 mixed interaction. It does not
identify path/state memory, noncommuting joins, Jordan structure, or a new
continuum field. The archive contains no ordered intermediate state.

A real context/morphism test still needs a semantics-matched
Gaussian-by-annulus `2 x 2` rectangle with all four cells in one random block,
and the tuple `(typed colour, ambient H1 before/after, J_local, R, ED/OD)` plus
one complete delete-one covariance. Existing Gaussian and annulus archives do
not provide synonymous rows, so no rectangle is fabricated here.
