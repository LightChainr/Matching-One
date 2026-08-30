# Minimal R-odd F3 charged source on Gaussian N13

All unweighted-null, charged-response and coefficientwise continuity gates pass.

## Charged basis

```text
q_A = T_01-T_10 = 1_axisX-1_axisY,
q_D = T_12-T_11 = 1_diagPlus-1_diagMinus.
```

Quarter-turn sends each charge to its negative. The projective action has C4 charge 2 and factors through C2; under full D4, q_A is B1 and q_D is B2.

## Explicit F3 defect

Reweight a rank-one state by `u^q_C`, with `u=omega`, `omega^3=1`. Then

```text
Z_C(u)=E[u^q_C],
O_C(u)=E[q_C u^q_C],
O_C(omega)=(omega-omega^2) W_C/2,
W_C=E[q_C^2].
```

Although the unweighted one-point is exactly zero, both charged responses are nonzero.

In the unit #337 H/A/D convention, the exact response matrix at `p_ref` is

```text
              unit A source   unit D source
H response       0               0
A response       0.159946977384               0
D response       0               0.0100607052647
```

| channel | W at p_ref | birth response | exit response | dW/dp |
|---|---:|---:|---:|---:|
| A / B1 | 0.319893954769 | 1.94280746874 | 1.84773531617 | 0.0950721525657 |
| D / B2 | 0.0201214105295 | 0.183555126861 | 0.214283068183 | -0.0307279413219 |

Every degree-12 current coefficient obeys `dW_C/dp=J_C,birth-J_C,exit`.

## Selection rule

An A source activates only A at linear order; a D source activates only D. H response and A-D cross-response remain exactly zero. The minimal production output is therefore `(W_A,J_A,birth,J_A,exit,W_D,J_D,birth,J_D,exit)` with one joint covariance block. Existing projective-line archives already contain it.

This constructs the missing explicit charged defect; it does not identify its continuum operator.
