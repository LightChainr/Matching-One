# Minimal lattice interface for the boundary CLE Q tangent

This is the Phase C bridge from the exact continuum tangent ODE to a finite
lattice measurement.  It deliberately uses plain connectivity indicators, so
the generic-(Q) projector problem is absent rather than hidden.

## 1. The three-component boundary vector is fixed by the source

Take four distinct boundary vertices (x_1<x_2<x_3<x_4) with free boundary
conditions.  The primary one-hot vector contains exactly the three link
patterns in Cai's Eq. (1.1):

\[
 \mathbf I=(I_{1234},I_{12|34},I_{14|23}).
\]

- `1234`: all four terminals belong to one FK cluster;
- `12|34`: (x_1,x_2) share one cluster and (x_3,x_4) a distinct cluster;
- `14|23`: (x_1,x_4) share one cluster and (x_2,x_3) a distinct cluster.

The separate `13|24` pairing is forbidden by planarity.  Configurations with
singletons or three-plus-one partitions are retained in the raw partition
label for normalization, but contribute zero to this three-vector.

Cai's Theorem 1.4 makes the useful branch identification exact:

\[
 U^{14|23}(\lambda)=C_1 V_{3h+1}(\lambda),\qquad
 U^{12|34}(\lambda)=C_1 V_{3h+1}(1-\lambda).
\]

Thus the high Frobenius branch from the parent artifact has a direct lattice
event, not a guessed projector.

## 2. Minimal sufficient statistics

On the square FK critical manifold

\[
 v=\sqrt Q,\qquad w_Q(A)=Q^{k(A)+b(A)/2},
\]

put (J=2k+b).  At (Q=1),

\[
 S_{\rm measure}=\frac{J-EJ}{2},\qquad
 \partial_Q E_Q[O]\big|_1=\operatorname{Cov}_1(O,J/2).
\]

Each synchronized batch and geometry stores only

```text
samples, sum_J, sum_J2,
for each of 1234,12|34,14|23:
    count, sum_J, sum_J2.
```

These integers reconstruct the ordinary probability, its measure tangent,
and the required delete-one covariance.  Recompute all ratios and covariances
inside every delete-one replicate.

## 3. Four derivative ledgers, never one catch-all score

For this observable the contributions are:

1. **Measure:** (\operatorname{Cov}(I_p,J/2)).
2. **Projector/bare field:** exactly zero, because the three terminal
   connectivity indicators have no explicit (Q) dependence.
3. **Boundary-field renormalization:** with lattice spacing
   (\delta=1/L), multiplying by (L^{4h(Q)}) contributes
   (4h'(1)\log L\,G), plus an unknown lambda-independent microscopic
   normalization tangent.
4. **Conformal prefactor:** Cai's Eq. (1.4) is

   \[
   G=K(x_i)^{2h}U(\lambda),
   \]

   so extracting (U) subtracts (2h'(1)\log K\,G).

Here (h'(1)=\sqrt3/(3\pi)).  Use normalized continuum coordinates for
(x_i) and keep (delta=1/L) separately; otherwise `log L` is counted both
in the field and in (K).

## 4. Four frozen cross ratios and the amplitude-free score

Freeze

\[
 \lambda\in\left\{\frac14,\frac13,\frac23,\frac34\right\}.
\]

A rational upper-half-plane representative is

\[
 (x_1,x_2,x_3,x_4)=
 \left(0,\frac{2\lambda}{1+\lambda},1,2\right).
\]

For each geometry and the `14|23` channel, form

\[
 z_i=rac{\partial_Q P_i}{P_i}
 +4h'\log L_i
 +\frac{\partial_Q^{\rm explicit}P_i}{P_i}
 -2h'\log K_i.
\]

The unknown derivative of (C_1(Q)) adds the same constant to every (z_i).
Subtract the value at (lambda=1/3).  This is the transparent version of
projecting the raw tangent off the ordinary-amplitude direction with the
covariance metric.

The exact-series parent oracle freezes the target

```text
lambda             1/4            1/3       2/3          3/4
d_Q log U anchored -0.1482999419   0          0.3249693721 0.3678625657
```

using 100 Frobenius coefficients; an independent 140-term evaluation changes
the values by less than `1e-12`.  The three non-anchor differences form the
residual (r).  Score

\[
 \chi^2=r^T\Sigma_r^+r,
\]

with degrees of freedom equal to the numerical covariance rank.  The
reflected `12|34` channel supplies a simultaneous crossing check on this
symmetric lambda grid; it is secondary, not extra independent evidence.

This is an inhomogeneous-ODE **solution residual**: the target is the uniquely
selected high-branch tangent solution of the differentiated ODE, modulo its
one amplitude gauge.  It avoids estimating third derivatives of noisy lattice
data by finite differences.

## 5. Tiny exact square-bond regression

The executable exhausts the 16 bond states of an open four-cycle with all four
vertices marked.  At (Q=1,v=1), it obtains

```text
pattern      probability   d_Q probability
1234         5/16          -37/256
12|34        1/16          -1/256
14|23        1/16          -1/256
total        7/16          -39/256
```

Every derivative computed from the integer sufficient statistics equals the
direct derivative of

\[
 \frac{\sum_A I_p(A)x^{J(A)}}{\sum_A x^{J(A)}}
\quad\text{at }x=1,\qquad \partial_Q=\tfrac12\partial_x.
\]

This tests the normalization and event semantics without claiming that a
four-edge graph is in the continuum regime.

## 6. No generic-Q sampler is needed

For plain connectivity probabilities, a (Q=1) FK stream plus (J/2) is the
exact tangent estimator.  If a color projector is added later, validate its
polynomial at integer (Q=2,3,4) by uniformly coloring FK clusters, then
differentiate the symbolic polynomial at (Q=1).  Do not estimate a derivative
by subtracting noisy simulations at different (Q).

## Claim layers

- Exact finite lattice: the (J/2) score, link-pattern classifier, derivative
  ledger, and four-cycle regression.
- Exact continuum input: Cai's branch identification and the parent branch's
  differentiated ODE/Frobenius tangent.
- Scaling hypothesis: square-FK boundary probabilities, after (L^{4h})
  renormalization, converge to the corresponding CLE connectivity vector.

## Reproduce

```bash
python3 scripts/p263_boundary_lattice_qscore.py \
  --output predictions/p263_boundary_lattice_qscore_20260829.json
python3 -m unittest discover -s tests -p 'test_p263_boundary_lattice_qscore.py'
```

## Primary source

- Gefei Cai, *Boundary four-point connectivities of conformal loop ensembles*,
  arXiv:2603.28161v2, <https://arxiv.org/abs/2603.28161>.

