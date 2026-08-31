# Two local sources: mixed response beyond the first-order Jacobian

## Result: own-source curvature hides an opposing normal response

The exact new64 tensor resolves own-source center curvature at both sizes
and orientations: `(-6.222±1.208,-5.985±1.148,-8.493±1.014,
-6.520±0.836)×10^-8`. The corresponding A(p_ref) curvatures are positive.
Mixed fs terms do not show a consistent cross-size direction. This
describes predominantly own-source bending in the specified coordinates,
not evidence for exact additivity.

The stronger new result comes from decomposing that bending. An exact
census projection removes the part aligned with the first two density
scores. The remaining own-center response is **positive** in all four
cases: `(4.116±0.949,3.233±0.669,3.300±0.954,3.977±0.682)×10^-8`.
The larger negative first-score tangent term and this positive normal term
partially cancel to produce the observed negative raw curvature.
[Source-normal construction and full result](p334-source-normal-curvature.md)
give an operational third kind of Euler-invisible perturbation, orthogonal
to the two original scores; this does not imply a third independent field.

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

## Why a rank-two Jacobian does not answer this question

An elementary finite-probability example separates the hypotheses exactly.
It is an abstract one-class example, not a constructed percolation lattice.
Let four equiprobable labels have `(L_f,L_s)` equal to the four binary pairs,
so `s_i=L_i-1/2` and the two-source tilted label law factorizes. Define

```
C_f = 1/2 + s_f/4 + gamma s_f s_s,
C_s = 1/2 + s_s/4,                 0 <= gamma <= 1/2.
```

Both centers lie inside[1/4,3/4]. One may give each a fixed small lifetime
and form valid ordered birth-clock pairs around these centers. At zero
source, the two-center Jacobian is `diag(1/16,1/16)` for **every** gamma;
the exact source Gram is also unchanged. Yet

```
H_fs C_f = gamma/16,   H_fs C_s = 0,
R_Cf = gamma tanh(1/4)^2,   R_Cs = 0
```

for the same +/-1/2 source rectangle. Thus identical source strength and
first-order local rank can coexist with either zero or nonzero mixed
response. The new coordinate probes how the future observable combines
the two inputs, even when the input law factorizes. It does not merely
repeat the first-order rank test.

## Provenance

The rectangle reader `c3a1b414` consumes saved20-batch finite-source vectors
only, with no finite-weight or trajectory recalculation. The resulting
[`score.json`](../results/p334-mixed-source-rectangle/score.json) retains
every corner, all physical/S/D readouts, cells, LOO and a factor on the
same deleted-batch sign convention used by the execution-team archive.
Original8 and new64 use the same prefix population; new64 only exists
on original00. They are not independent population replications.

The second-score tensor is delivered at
[`c48fa360`](https://github.com/LightChainr/Matching-One/blob/c48fa360a37a9887ef32ff6d3ce947c4e4601b53/notes/p334-mixed-source-curvature.md).
The exact policy/Hessian/rectangle identities and symmetry analysis are
at [`a6e7141a`](https://github.com/LightChainr/Matching-One/blob/a6e7141ac4a7a0b34fca26373963681aa8534de6/notes/p334-fixed-source-mixed-curvature.md).
