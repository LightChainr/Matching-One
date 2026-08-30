# P250: the natural physical site toggles cannot carry order

Date: 2026-08-30. Status: **exact obstruction; no pilot runner added**.

The requested physical embedding of the P333 detach/join rectangle does not exist in the current P215/P225/P334 runners.  This is now an exact operator statement rather than a missing-interface impression.

## Runner audit

The three relevant lineages all use site overwrite maps:

- P215 forces one fixed site absent and then present on a frozen field.  Its parity companion evaluates the same insertion in black NN and complement-reversed white matching connectivity.
- P225 forces the pivotal root absent/present and compares the two wrapping outcomes.  Multiradius changes the landing measurement, not the mutation.
- P334 adds the next absent vertex in a fixed permutation.  Its `K1/K2` and first-line records retain path timing, but it does not define a delete/join pair or emit a common four-corner rectangle.

There is no existing state-dependent connector selection whose microscopic support can change after the first operation.

## Exact obstruction

Let `F_v^0` force site `v` vacant and `F_v^1` force it occupied.  Write

```text
D = F_d^0,
J = F_j^1,
```

and apply words left-to-right.

If `d != j`, the maps act on disjoint Boolean coordinates, hence

```text
DJ = JD
```

as complete site fields.  Consequently every deterministic response, including a typed projective leg, obeys

```text
R_minus = L_DJ - L_JD = 0.
```

This remains true even though the symmetric mixed difference

```text
R_plus = L_D + L_J - L_DJ - L_JD
```

can be nonzero.

If `d=j`, overwrite absorption gives

```text
DJ=J,  JD=D,
```

so `R_plus=0` for every response.  Moreover one base site cannot simultaneously satisfy the proposed physical preconditions “occupied site to delete” and “vacant connector to add.”

Thus no fixed-site choice makes both the connected P333 channel and the P250 order channel nondegenerate.

## Smallest honest square-torus certificate

The script exhausts all 64,512 admissible `(field,root,d,j)` rectangles on the `L=3` square torus, with three distinct marks, `d` occupied and `j` vacant.  `L=3` is used because `L=2` has collided nearest-neighbour images.

For the production-type response

```text
black-NN rank-one root membership
minus white-matching rank-one root membership,
```

the exact histograms are

```text
R_plus:  -2:4212, -1:15696, 0:25956, 1:16056, 2:2592
R_minus:  0:64512.
```

A minimal transparent `R_plus=1` witness uses the vertical sites

```text
root=(0,0), d=(0,1), j=(0,2), base occupied={(0,0),(0,1)}.
```

Adding `j` closes the vertical rank-one loop, whereas deleting `d` first prevents that closure.  The two final fields after `DJ` and `JD` are nevertheless identical, so `R_minus=0`.

The same enumeration verifies the exact typed involution

```text
L_NN(omega,r) = -L_matching(C omega,r),
C D_v = J_v C,
C J_v = D_v C.
```

So simultaneous black-NN/white-matching bookkeeping is not the problem; loss of ordered state is.

## Decision and minimum escape

No pilot interface is added.  Renaming fixed overwrites as formal P333 detach/join would create a guaranteed-zero order channel.

A future physical rectangle must first declare a genuinely state-dependent cut/connector morphism.  At minimum it must save:

- the connector/cut support selected separately after each first move;
- black-NN and white-matching component IDs, ambient ranks, and primitive lines in `S_D,S_J`;
- `L_D,L_J,L_DJ,L_JD` on the same base field and marked triple.

The support itself is essential: without it, a nonzero response cannot certify that both branches applied the intended morphism.  This note does not rule out such a declared adaptive operation; it rules out the natural fixed-site embedding already present in the repository.

## Reproduction

```sh
python3 scripts/p250_physical_toggle_obstruction.py
python3 tests/test_p250_physical_toggle_obstruction.py
```

Output: `results/exact-p250-physical-toggle-obstruction/latest.json`.
