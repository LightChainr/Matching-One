# P537 cyclic bridge-order thermal rotation

## New exact result

The two surviving six-block bridge-order modes in the provisional N25 landing
contract do not merely have different amplitudes.  Thermal differentiation
rotates them differently in the root-conditioned P4 plane.

Write

```text
F_s = jY_s - (Y_p/M_p) jM_s
H_s = T_{p,s}/M_p = (dF_s/dp)/M_p
```

for source mode `s`.  The exact interval matrix with columns
`(clean_same, clean_reversed)` is numerically

```text
[[ +4.643561586827378e-08, -4.263577227418322e-06 ],
 [ +7.600979820643663e-07, +5.046234984612301e-06 ]].
```

Its determinant is

```text
+3.475061476262754e-12
```

and its exact rational interval excludes zero.  Equivalently, the normalized
thermal slopes have different signs:

```text
H_same/F_same         = +16.3688575644...
H_reversed/F_reversed =  -1.1835683313...
```

Thus the root-conditioned thermal jet resolves a genuine two-dimensional
cyclic-order response.  A scalar “clean bridge amplitude” cannot absorb both
orders: the `clean_same` mode starts with a small positive fixed-M response
and grows in the positive thermal direction, while `clean_reversed` starts
negative but has a positive thermal derivative.

## Scientific use

This identifies the cyclic order of the two bridges as a dynamical coordinate,
not just a combinatorial subdivision of one landing amplitude.  Any reduced
ordinary-channel model that merges the two orders before taking the thermal
jet loses an exact rank-two structure.  The natural next question is whether
the same two directions survive the complete site-flip landing tensor, where
the pivot, rank transition, source-absent Schur term, and off-port branch and
separator identities are retained.

## Boundary

This is an exact post-processing theorem about the same N25 populations used
by `p537-landing-matrix-preflight-20260901`.  It introduces no sampling error,
but it is not independent evidence.  More importantly, its Bell-port
`clean_two_bridge_six_block` event is still provisional: it has no explicit
thermal pivot and cannot certify the repository's intended
ordinary-site-flip/no-extra-branch event.

Run:

```bash
python3 experiments/p537-cyclic-bridge-jet-20260901/score_cross_mode_jet.py
```
