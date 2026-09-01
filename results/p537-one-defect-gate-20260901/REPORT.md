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

## Distance-two typed-carrier supplement

The frozen row-major configuration `(x,y,z)=(0,6,2)` with off-z occupied
cells `{1,3,4,5,7,9,10,12,15,16,17}` gives a second exact edge:

```text
internal centers             x=0, y=4, z=3
pairwise NN distances        2, 2, 2
rank index                   0 -> 1
Bell                         274568 -> 8256
joint C                      21990249529872 -> 3298535014656
joint terminal incidence     2 -> 1
g16                          8 -> 0
Delta a                      -1/50
```

Its globally frozen pooled-root C4-orbit weight is

```text
source midpoint       -1.0121115955209059e-10
root counterterm      +5.341390620686106e-12
full                  -9.586976893140449e-11
```

All three signs exclude zero by exact rational intervals.  The source touches
the global black landing component and the common typed carrier changes even
though no pair of marked centers is within NN distance one.  Hence the metric
split `d_NN <= 1` cannot define the full contact/collision contribution.  The
correct finite classification is joint-incidence/typed carrier.  Pairwise
distance two on N5 is not an annular certificate, so the separated sector and
any asymptotic consequence remain open.
