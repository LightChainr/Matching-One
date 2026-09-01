# P537 one-defect diagonal-edge result

The frozen axis-N25 scan stopped after 12,568 backgrounds at the first
literal site flip that changes rank, the source Bell partition, and the
actual source value:

```text
transition       axis-N25:x0:zE:y(-1,-1):eta12567:0to1
rank index       0 -> 1
Bell             9240712 -> 6848576
joint C          23090870354448 -> 92359816642816
g16              4 -> 0
Delta a          -1/100
```

At the globally frozen two-geometry pooled root, the restored C4-orbit and
geometry-pool-half edge weight is

```text
source midpoint       -1.0888815582478189e-11
root counterterm      +2.5901918540035547e-12
full                  -8.298623728474635e-12
```

Every displayed sign has a rational interval certificate in `result.json`.
The beta-free source term already excludes zero, so the decision is robust to
the concern that a Bell-cell counterterm allocation created the witness.

This is a contact/collision edge: `arm_mask=3`, the source touches the global
black landing component, and the alternating-collar identity does not apply.
It falsifies an automatic two-independent-defect/six-arm gain for the full
physical graph.  It does not establish failure in the ordinary separated
sector or an asymptotic lower bound.

Reproduce all committed outputs with:

```bash
python3 experiments/p537-one-defect-gate-20260901/verify.py
```
