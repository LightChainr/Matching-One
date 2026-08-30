# P333/P321/P370 minimal multimark lower-bound gate

Neither lower bound from 081a5ed is attained by the frozen falling-factorial marks. Width three with two marks has ladder 10 -> 8 -> empty; width four with three marks has 17 -> 15 -> empty. In both cases the Gram restriction has coefficient rank zero and augmented rank one, so every surviving affine modulus is invisible to the obstruction. The family-specific bounds rise to three total marks at width three and five at width four. The width-four bound exceeds its four nonzero falling-factorial responses, exhausting that scalar family without testing another case.

| width | marks | dim W | rank G0 | dim radical | mark velocity rank | affine jet | + endpoint/radical | + Gram | + source | decision |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 3 | 2 | 7 | 5 | 2 | 1 | 10 | 8 | empty | empty | `lower_bound_not_attained` |
| 4 | 3 | 17 | 7 | 10 | 2 | 17 | 15 | empty | empty | `lower_bound_not_attained` |

## Exact interpretation

- Width 3 with exactly 2 marks: `lower_bound_not_attained`; canonical Gram-skew rank 2.
  First empty restriction `endpoint_radical_normalized -> gram_self_adjoint` with `y^T C=0 but y^T b=1`. The family-specific lower bound rises to 3 total marks; available nonzero falling-factorial marks=3, family exhausted=False. The larger case was not tested.
- Width 4 with exactly 3 marks: `lower_bound_not_attained`; canonical Gram-skew rank 4.
  First empty restriction `endpoint_radical_normalized -> gram_self_adjoint` with `y^T C=0 but y^T b=1`. The family-specific lower bound rises to 5 total marks; available nonzero falling-factorial marks=4, family exhausted=True. The larger case was not tested.

## Boundary

- Exact rational Q=1 first-jet algebra only for the two frozen cases.
- Marks retain only falling-factorial endpoint responses, not rooted-cluster geometry or arbitrary matrix degrees of freedom.
- No continuum LCFT, physical transfer matrix or formal-semigroup K identification is made.
- If the bound fails, only a new lower bound is reported; no additional mark count is tested.
