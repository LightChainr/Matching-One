# Phase E: path-ordered filtration witness

**Mechanism: the two joins have the same endpoint, but they can expose a final
rank-two homology event at different stopping times.**  This is the first
Issue #200 observable that is genuinely ordered without contradicting the CRT
join theorem.  It uses an intermediate sigma-algebra; it is not a function of
the final join alone.

The oracle exhausts the same honest Gaussian lift as Phase C,

```text
N1 -> N2/N5 -> N10,       (1+i)(2-i)=3+i,
P(3+i) = [[3,-1],[1,3]],  column-HNF [[10,3],[0,1]].
```

All `2^10=1024` masks are evaluated with black-NN and white-matching layers,
real quotient keys, and exact displacement-aware union-find.  The intermediate
HNFs are `[[2,1],[0,1]]` and `[[5,3],[0,1]]`.  As in Phase D, ambient H1 after
an artificial fibre identification is labelled by the frozen raw
representative-displacement convention.

## Ordered rank filtration

For one typed colour let `h0,h2,h5,h25` be its ambient-H1 ranks at the four
join corners, and set

```text
A = 1{h25=2}.
```

Compare the two finite filtrations

```text
F(2,5): sigma(h0) < sigma(h0,h2) < sigma(h0,h2,h5,h25)
F(5,2): sigma(h0) < sigma(h0,h5) < sigma(h0,h2,h5,h25).
```

Their endpoint sigma-algebras are identical.  Conditional on `A=1`, define
the first stage `tau` in `{0,1,2}` where rank two is visible.  The primary toy
path mark is

```text
C_act = tau(5,2)-tau(2,5)
      = A [1{h2=2}-1{h5=2}].
```

This mark is odd under path reversal `R2 <-> R5`.  It cannot be reconstructed
from `(h0,h25)`: it explicitly asks what the first intermediate rank was.

## Exact nonzero witness and typing

Balanced mask `62`, occupied sites `[1,2,3,4,5]`, has

```text
black NN:       (h0,h2,h5,h25)=(1,1,1,1), C_act= 0
white matching: (h0,h2,h5,h25)=(1,1,2,2), C_act=-1.
```

Thus its typed-complement rows are

```text
C_even=(C_B+C_W)/2=-1/2,
C_odd =(C_B-C_W)/2=+1/2.
```

Over all masks at `p=1/2`, the exact means are

```text
E C_B    = -75/1024,       nonzero on 75 masks
E C_W    = -21/512,        nonzero on 50 masks
E C_even = -117/2048,
E C_odd  = -33/2048.
```

The `even/odd` labels refer to the exact typed layer swap
`black-NN <-> white-matching`.  `C_act` is path-orientation odd.  This one-fibre
N10 oracle does **not** determine the geometric N650 `S/D` parity; that needs
both real N650 geometries or a second exact finite geometry.

## Doob and information checks

For the same final event, the oracle computes

```text
C_Q = sum (Delta E[A|F_k(2,5)])^2
    - sum (Delta E[A|F_k(5,2)])^2.
```

It is configurationwise nonzero on `669` black masks and `382` white masks.
Masks `141` and `173` have the same black endpoint ranks `(h0,h25)=(0,2)` but
different `C_Q`, proving that the result did not collapse to a final-join
function.  Its uniform mean is exactly zero in both layers, as required by
martingale isometry; the zero mean is a theorem, not evidence that the paths
carry the same intermediate information.

The conditional-variance diagnostic is nonzero in mean:

```text
E[Var(A|F1(5,2))-Var(A|F1(2,5))]
  = -110932591235/1819605946368   black NN
  = -6026549/190440448            white matching.
```

Negative sign means that, in this exact convention, the `R5` first-stage rank
is more informative about the final rank-two event than the `R2` first-stage
rank.

## Rc scalar-rank no-go and minimum retained state

Repeating the Doob construction for the Phase D target
`1{R_nonlocal != 0}` with the natural scalar summaries

```text
(r0), (r0,r2)/(r0,r5), (r0,r2,r5,r25,J_local)
```

gives `C_Q=0` on every mask and both typed layers.  So the enormous Phase D
`Rc` rejection is a symmetric nonlocal interaction and still cannot be called
chronological memory.

The smallest future acquisition for the nonzero H1 path mark is to retain
per-colour `h0,h2,h5,h25` before batch aggregation.  To make `Rc` itself
chronological requires strictly richer data: a dynamic edge order or marked
cluster lineage at the intermediate join.  Scalar ranks are exactly
insufficient in this oracle.

The frozen 100M Phase D job remains stopped.  Phase E is exact/toy evidence
only and makes no Jordan or continuum-memory claim.

Reproduce:

```bash
python3 scripts/p200_path_filtration_oracle.py \
  --output results/exact-cover-character-oracles/p200_path_filtration.json
python3 -m unittest discover -s tests \
  -p 'test_p200_path_filtration_oracle.py'
```
