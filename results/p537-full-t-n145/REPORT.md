# Held-out N145 complete thermal transport

Status: completed at the frozen sample count; no top-up.

The independent N145 block used `200,000,000` common-field samples, 16
deterministic shards, 100 source batches and seed `20260901545`.  It was scored
jointly with the independent 100-batch P50 N145 full-curve baseline.  The
complete canonical source and Schur/root projection give

```text
p_root = 0.592741991414186 +/- 0.000003168661463
T_t    = -1.16386223543e-6 +/- 1.85491364273e-7
J145   = -0.000694364434794 +/- 0.000110663607588
95% CI = [-0.000911265105666, -0.000477463763921]
```

## Frozen decision

The frozen boundary was `J=-0.0004931`.  The upper confidence endpoint lies
only `1.5636e-5`, or `0.1413` standard errors, above that boundary.  The
interval therefore narrowly straddles it and the literal frozen outcome is

```text
N145_FULL_T_UNRESOLVED
```

No sample extension is allowed or performed.

## Scientific discrimination retained by the result

The unresolved stop rule does not make the observation neutral.  Using the
exact N25 complete `J=-0.00551943142484`, the held-out N145 center gives the
two-size power

```text
q_J(25,145) = 1.17930 +/- 0.09066.
```

The fixed `q_J=5/4` law predicts `J145=-0.000613210171405`, only `0.733`
N145 standard errors from the observation.  The fixed faster `q_J=7/4` law
predicts `-0.000254621664678`, `3.974` standard errors away.  This comparison
was not substituted for the preregistered interval gate, but it does eliminate
the faster-decay model as a useful description of the N25-to-N145 full
transport while leaving the four-packet power compatible.

This is a full-`T` result and is invariant under the common thermal coordinate
shift.  The separate contact-stage gauge question is handled by the horizontal
sidecar, not by altering this held-out score.
