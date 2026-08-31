# Two local sources: mixed response beyond the first-order Jacobian

## The next mechanism question

Two locally independent first-order birth-center response directions have
already been measured, and exact contact features organize much of their
prefix variation. This readout asks a different question: are the two knobs
additive under their **specified finite source policy**, or is there a
measurable mixed response? No predictive model or conditional-shape test is
repeated here.

The two-source family holds every joint-safe degree class mass fixed:

```
q(t_f,t_s;u|Z) = pi_a exp[pi_a(t_f L_f(u)+t_s L_s(u))]
                 / sum_{v in a} exp[pi_a(t_f L_f(v)+t_s L_s(v))].
```

Outside safe classes q=1/d. Source marks are the same R0-only loop counts
as before. Both geometries' immediate rank/Euler joint distribution is
unchanged for every finite pair(t_f,t_s).

## A rectangle already present in the saved data

The previous `g_plus=(L_f+L_s)/2` and `g_minus=(L_f-L_s)/2` runs at t=±1
are exactly the four physical source corners(±1/2,±1/2). Their mixed contrast

```
R_F = Delta_plus(+1)+Delta_plus(-1)
      -Delta_minus(+1)-Delta_minus(-1)
```

is the integral of `partial_f partial_s E_q F` over this unit-area square.
It equals an average mixed derivative across the box, not necessarily the
derivative at zero. For any fixed-coordinate additive law
`E_q F=A(t_f)+B(t_s)`, it is zero. Every original prefix outside00 also
has exactly zero mixed contrast because at least one source is identically
absent. No source or geometric-exchange symmetry is assumed within00.

Saved finite responses from `8ad30617` are changes from a common baseline,
which cancels in the rectangle. C responses are minus one half of the
integrated A responses; W responses are minus integrated E responses.
The original full20000-prefix denominator and20 paired batches remain.

The initial original8 rectangle is weak: most physical mixed responses are
within approximately2 SE, with different signs between geometries or sizes.
For example mixed C is `(-2.970±2.780,3.878±2.519)e-8` at N325 and
`(2.696±2.805,4.651±2.273)e-8` at N425, in first/second order. It neither
establishes additivity nor a consistent nonlinear coupling. The exact
second-score readout from the already collected new64 is the more precise
next coordinate, without allocating new trajectories.

## Exact zero-source second derivative

With `s_i=pi_a(L_i-mu_ai)`, define

```
t_ij = s_i s_j - pi_a^2 Cov_a(L_i,L_j).
H_ij F = E_uniform[t_ij F].
```

The subtraction matters: multiplying two first scores alone changes class
mass and is not the second derivative of the declared normalized policy.
Let n_a be class count, S_i the summed marks and Q_ij the summed product.
The exact integer numerator of t_ij, with denominator d², is
`(n_a L_i-S_i)(n_a L_j-S_j)-(n_a Q_ij-S_i S_j)`.
Its class-weighted sum vanishes exactly. For independent next labels U,V,

```
H_ij F = E[(t_ij(U)-t_ij(V))*(F_U-F_V)/2].
```

F_U averages the two existing independent suffixes at label U. Retaining
the ff/fs/ss tensor also separates mixed coupling from each source's own
curvature. All four observables A(p_ref),E(p_ref),C,W remain paired.

The natural source coordinate choice is part of the physical protocol.
A nonlinear reparameterization can change raw Hessian components; they
are not continuum-field identities or coordinate-free manifold curvature.
Mixed derivatives commute here. A nonzero entry does not demonstrate
temporal memory or a noncommuting order of operations.

## Provenance

The rectangle reader `c3a1b414` consumes saved20-batch finite-source vectors
only, with no finite-weight or trajectory recalculation. The resulting
[`score.json`](../results/p334-mixed-source-rectangle/score.json) retains
every corner, all physical/S/D readouts, cells, LOO and a factor on the
same deleted-batch sign convention used by the execution-team archive.
Original8 and new64 use the same prefix population; new64 only exists
on original00. They are not independent population replications.
