# Minimal explicit R-odd F3 charged source

Status: exact N13 state and directed-boundary response on the Gaussian quotient
`3+2i`.  This follows the unweighted A/D null theorem in `770fe27` by inserting
the smallest source that transforms in the missing quarter-turn-odd sector.

## Representation and source

Order the four F3 projective lines as

```text
(axis x, axis y, diagonal +, diagonal -).
```

The physical quarter-turn swaps the two entries in each pair.  Hence the
projective permutation representation splits into two even orbit sums and the
two charged directions

```text
q_A = 1_axisX - 1_axisY = T_01-T_10,
q_D = 1_diagPlus - 1_diagMinus = T_12-T_11.
```

Both have quarter-turn eigenvalue `-1`.  This is C4 charge 2, not a `+/-i`
doublet: `R^2=-I` acts trivially on projective lines, so the action factors
through C2.  With the declared axis reflection, `q_A` is the D4 `B1` block and
`q_D` the `B2` block.

For either charge `q_C in {-1,0,1}`, insert the explicit defect

```text
state weight = exp(s_C q_C)
```

or its genuinely F3-valued phase specialization

```text
u=omega, omega^3=1, state weight=u^q_C.
```

Define

```text
Z_C(u)=E[u^q_C],
O_C(u)=E[q_C u^q_C]=u dZ_C/du,
W_C=E[q_C^2].
```

Exact pair symmetry gives

```text
Z_C(u)=1-W_C+(u+u^-1)W_C/2,
O_C(u)=(u-u^-1)W_C/2.
```

Thus the unweighted one-point still vanishes, while the charged one-point at
`u=omega` and the zero-source susceptibility are nonzero.  This is the minimal
explicit way to activate an R-odd sector without breaking the quotient graph
or inventing a new path label.

## Selection matrix

An A source activates A and a D source activates D.  The H response is zero by
R parity.  A-D mixing is also exactly zero because `q_A q_D=0` statewise.  In
the normalized #337 coordinates `A=q_A/sqrt(2)`, `D=q_D/sqrt(2)`, the response
matrix at `p_ref` is

```text
              unit A source   unit D source
H response       0                 0
A response       0.159946977384    0
D response       0                 0.0100607052647
```

The corresponding integer-charge susceptibilities and currents are

| sector | `W_C` | birth | exit | `dW_C/dp` |
|---|---:|---:|---:|---:|
| A / B1 | 0.319893954769 | 1.94280746874 | 1.84773531617 | +0.0950721525657 |
| D / B2 | 0.0201214105295 | 0.183555126861 | 0.214283068183 | -0.0307279413219 |

The D response is about sixteen times smaller than A at this reference point,
but is exactly nonzero.  Its negative slope says the diagonal charged weight is
already leaving faster than it is born there; this is a finite-volume response
fingerprint, not yet a scaling claim.

## State, birth and exit oracle

For every microcanonical size `k`, the exact subset oracle records the two
line counts in the selected orbit.  Their sum divided by `C(13,k)` is the
Bernstein coefficient of `W_C`.  At every lower size `k=0,...,12`, the sums of
the corresponding directed birth and exit edges divided by `C(12,k)` are the
charged currents.  The certificate verifies coefficient by coefficient that

```text
dW_C/dp = J_C,birth - J_C,exit
```

for both A and D.  The same row also contains the Laurent coefficients of
`Z_C(u)` and `O_C(u)`, so state, source and sink response are independently
recoverable rather than inferred only at one value of `p`.

## Minimal next production quantity

The existing projective-line archive is already sufficient.  The production
tuple should be

```text
(W_A, J_A,birth, J_A,exit, W_D, J_D,birth, J_D,exit)
```

from the same random block with its complete cross-channel covariance.  No new
`iota`, path-memory, or marked-site field is required.  If one complex output
is preferred, publish `O_A(omega)` and `O_D(omega)`; each is a fixed phase times
`W_C`, so no phase fit is allowed.

## Boundary

This is an exact finite charged-defect construction and a complete linear
selection rule for the two R-odd projective multiplicities.  It does not prove
that either response survives a continuum limit, assign a CFT primary, or turn
the former unweighted N65 diagonal-odd fluctuation into evidence.

Reproduce with:

```bash
python3 scripts/p337_n13_r_odd_charged_source.py \
  --json results/p337-n13-r-odd-charged-source/latest.json \
  --markdown results/p337-n13-r-odd-charged-source/latest.md

python3 -m unittest discover -s tests \
  -p 'test_p337_n13_r_odd_charged_source.py'
```
