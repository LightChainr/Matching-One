# P337 thermal/pivotal gate audit: the covariance jet has two channels

## Decision

The exact-​`p_c` summability theorem for the raw canonical pair kernel does
not by itself control the thermal derivative that enters original `U`.  The
missing object is the derivative of a covariance, and it contains two terms.
Controlling only `partial_p E[g_xy]` drops the topology-readout pivotal term
and can even give the wrong sign on finite controls.

This audit therefore retires the shortcut “prove raw-kernel derivative
summability, then declare the original-`U` gate closed.”  It does not prove
that the complete thermal response grows or vanishes.  Issue #537 is the one
canonical theory task for that remaining question; no new Monte Carlo block is
authorized.

## Exact finite identity

For occupation functions `O` and `a` without explicit `p` dependence, write

```text
D_z f     = f(z=1)-f(z=0),
f_z^mid   = [f(z=1)+f(z=0)]/2.
```

Under the Bernoulli law on the other sites,

```text
d_p Cov_p(O,a)
 = sum_z E_-z[(O_z^mid-E_p O) D_z a]
 + sum_z E_-z[(a_z^mid-E_p a) D_z O].
```

The first term is the canonical-kernel reconnection channel.  The second is
the original rank/readout pivotal channel.  They are an algebraic
decomposition of one derivative, not two independent evidence blocks or
causal percentages.  Direct probability differentiation, the midpoint
identity and the Bernoulli score agree exactly in every finite control below.

For the spatially averaged canonical source `a=N^-2 sum_(x!=y) g_xy`, the
endpoint contribution to the raw mean derivative is also exact:

```text
sum_(z=x,y) E D_z g_xy = -2 E[g_xy]/(1-p).
```

This explicit vacancy-dilution term is already controlled by raw-kernel
summability and is not a long-range pivotal enhancement.

## What one ordinary site flip can do

For `z` distinct from the two marked vacancies, fill `z` after forming the
physical-edge equality components with `z` vacant.  Filling `z` merges at most
the four components touching its incident edges.  On the marked eight ports,

```text
D_z g_xy = g(join_T pi_xy)-g(pi_xy),   |T|<=4.
```

The complete Bell-8 audit contains 4,140 partitions and 64,954 allowed joins.
It finds 29,970 nonzero differences: 9,952 positive and 20,018 negative, with
the exact abstract envelope `|16 D_z g|<=68`.  A nonzero difference never
starts from zero shared terminal components and always has at least two shared
components on one side of the flip.  Shared-component count alone is still
insufficient: 312 actual N13 flip/pair records preserve that count while a
kernel sign changes.

The full N9/N10/N13 torus enumerations cover 9,728 configurations and, for the
N13 full pair/flip check, 1,198,080 state-edge-pair records.  Lifted topology,
digital-dual rank and the port-join update agree exactly.  At `p=3/5`, `O=E`,
the two covariance-derivative terms are

| N | kernel reconnection | readout pivotal | total |
|---:|---:|---:|---:|
| 9 | +0.0076368394584064 | -0.0248721173905408 | -0.0172352779321344 |
| 10 | +0.0117678699970560 | -0.0255705979992146 | -0.0138027280021586 |
| 13 | +0.0061685367582240 | -0.0166271790220691 | -0.0104586422638452 |

Omitting the readout term reverses the sign in all three controls.  These are
finite semantic controls, not critical-scale evidence.

## The actual asymptotic gate

At the pooled root, retain the repository notation

```text
M=mean_g <q>,       Y=P4 <E>,       D=M_p,
R=Y_p/D,            jM=mean_g Cov(q,a),
jY=P4 Cov(E,a),     A_N=N^(13/8)/2.
```

The quantity that must be controlled is

```text
T_N = jY_p - R*jM_p - R_p*jM,
J_N = A_N*T_N/D.
```

To prove this source disappears from original `U`, one needs
`T_N=o(D/A_N)`.  Raw summability, thermal summability without a rate, or a
positive finite `D` does not imply this.  Under the conditional scale
`D~N^(3/8)`, the required bound is `T_N=o(N^-5/4)`.

A naive absolute three-position estimate is not enough with the currently
available rigorous square-site arm exponent.  If comparable-scale `x,y,z`
require three alternating four-arm events, one dyadic shell has the form

```text
R^4 * pi_4(R)^3.
```

Absolute summability by this route needs an exponent strictly larger than
`4/3`; the rigorous input used for the raw pair theorem gives only some
`alpha_4>1`.  Even inserting the unproved square-lattice value `5/4` would
leave `R^(1/4)`.  The remaining work must therefore find stronger geometric
support, a signed cancellation/landing identity, or a different normalization
estimate.  It must also transport exact-`p_c` control to the pooled root.

## Reproduction and scope

The audit is in
`experiments/p337-thermal-gate-audit-20260901/thermal_gate.py`; it uses only the
Python standard library and refuses geometries above N13.  It is a bounded
oracle, not a production engine.

```bash
python3 -m unittest discover \
  -s experiments/p337-thermal-gate-audit-20260901 \
  -p 'test_thermal_gate.py' -v
python3 experiments/p337-thermal-gate-audit-20260901/thermal_gate.py \
  --join-audit --out /tmp/p337-joins.json
python3 experiments/p337-thermal-gate-audit-20260901/thermal_gate.py \
  --geometry 3 2 --out /tmp/p337-n13.json
```

The 12 focused tests pass.  No new random configuration, distance, exponent
fit, completion scan, cloud job or GPU task was used.  A full N25
kernel/readout split can remain P2 reproducibility support, but reproducing the
already known `J2` does not close this asymptotic gate.
