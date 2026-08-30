# P337/P334: direct rank-two births have a theta/figure-eight carrier

Date: 2026-08-30. Status: **deterministic finite-torus lemma plus exact bounded census**.

This note closes the first geometric gap behind the direct-birth mass `D_N=P(K1=K2)`.  It does not import a six-arm exponent and does not identify a CFT field.

## Result

Let `T=Z^2/Lambda` carry nearest-neighbour square-site connectivity.  Let `S` be occupied, let `v` be vacant, and assume the ambient image of `H1(S)` in `H1(T)=Lambda` is zero.  For each occupied component `C` of `S`, choose one component `C_tilde` of its inverse image in `Z^2`.  Every occupied neighbour edge `e=(v,u)` has a lift adjacent to a fixed lift `v_tilde`; that lifted endpoint lies in a unique translate

```text
lambda_e + C_tilde,   lambda_e in Lambda.
```

Then the ambient homology created by adding `v` is exactly

```text
L_v(S) = sum_C span_Z { lambda_e-lambda_f : e,f meet the same old component C }.
```

Consequently a direct jump `0 -> 2` on the degree-four square lattice has only two carrier types:

1. **one-carrier theta:** one old component has deck addresses of affine rank two.  Three of its addresses are affinely independent;
2. **two-carrier figure-eight:** no old component has affine rank two.  Degree four then forces two old components, each met twice, and their two nonzero address differences are independent.

This is the exact replacement for the informal claim that every direct birth is simply “a six-arm site”.  The first type has at least three long occupied carrier arms.  The second has at least four and is an eight-arm, not an exactly-six-arm, local topology.

### Proof of the cover-address identity

Since the target deck group is abelian, zero ambient `H1` image also makes the image of every old component's fundamental group trivial.  Its inverse image is therefore a disjoint union of copies mapped homeomorphically onto that component.

Every new cycle alternates star edges at `v` with old paths inside occupied components.  An old path joining incidences `e,f` lifts between the adjacent endpoint at address `lambda_e` and a translate of the endpoint at address `lambda_f`; its deck displacement is `lambda_e-lambda_f`.  This proves containment in the displayed subgroup.  Conversely, the old path followed by the two star edges is a new cycle with exactly that displacement, proving equality.

If the subgroup has rank two, either one address set already has affine rank two, or two rank-one address sets supply independent directions.  In the latter case each set needs two incidences, so degree four leaves exactly the `2+2` figure-eight partition.  This proves the dichotomy.

## What is true about arms

Let

```text
s_infinity = min_{0 != lambda in Lambda} ||lambda||_infinity.
```

Fix an integer radius `R` for which the one-step enlarged square box is embedded, for example `2(R+1) < s_infinity`.  Distinct carrier addresses give distinct occupied components in the universal cover.  Each such lifted component contains both an endpoint one step from `v_tilde` and another endpoint within one step of a nonzero deck translate of `v_tilde`.  It must therefore cross the annulus from the star of `v_tilde` to radius `R`.

Thus the theta case supplies three disjoint occupied NN crossings and the figure-eight case supplies four.  In the planar embedded annulus, consecutive distinct occupied crossing components are separated by vacant crossings on the square-site matching lattice (NN+NNN).  Applying the elementary site-matching separation in every cyclic gap gives respectively

```text
theta:       3 occupied + 3 vacant alternating arms,
figure-eight: 4 occupied + 4 vacant alternating arms.
```

The extra one-step margin in the embedding condition accommodates the diagonal matching edges.  This is a microscopic-to-injectivity-scale deterministic implication, not an asymptotic probability equivalence.

## Minimal obstruction to an “exactly six” dictionary

On the `3 x 3` torus (`a=3,b=0,N=9`), use the repository's BFS vertex labels.  Take

```text
S mask = 30, occupied vertices = [1,2,3,4], birth site v = 0.
```

Before the birth, the horizontal pair and vertical pair are two rank-zero components.  Their deck-address differences are `(-3,0)` and `(0,-3)`.  Adding the origin closes two independent loops at once.  This is a figure-eight carrier, so the naive one-to-one identification with an exactly-six-arm theta event is false.

It is the first figure-eight in the bounded Gaussian census used here: norm 5 has no direct birth; norm 8 has 40 direct births, all theta; norm 9 has 36 theta and 9 figure-eight births.

## Exact census

The certificate exhausts every direct subset edge on the requested small quotients.

| generator | N | direct edges | theta | figure-eight |
|---|---:|---:|---:|---:|
| `2+i` | 5 | 0 | 0 | 0 |
| `2+2i` | 8 | 40 | 40 | 0 |
| `3` | 9 | 45 | 36 | 9 |
| `3+i` | 10 | 80 | 80 | 0 |
| `3+2i` | 13 | 793 | 793 | 0 |
| `4` | 16 | 4,624 | 4,288 | 336 |
| `4+i` | 17 | 8,823 | 8,704 | 119 |

The N=9,10,13,16,17 totals reproduce the existing direct-edge table.  Every edge is classified by exact integer lift coordinates; no geometric tolerance is used.

## Consequence for the `D_N` program

The deterministic lemma gives only the easy inclusion

```text
direct 0->2 birth  subset  polychromatic six-arm event to injectivity scale.
```

It therefore supports an **upper comparison** with a six-arm probability.  A matching lower comparison needs a genuinely new nondegeneracy/gluing statement: conditional on the six arms, the three occupied arms must close into two independent deck directions with probability bounded away from zero, uniformly in scale and in the near-critical window.  One must also control the figure-eight share, window tails, and any lattice-universality transfer.  Without those ingredients, `D_N ~ N^(-5/6)` remains a conditional heuristic rather than a theorem.

The converse deterministic implication is false: six alternating arms can live inside a contractible region and need not create either torus generator.  Hence direct births form a globally typed subset of the ordinary six-arm event.

## Reproduction

```sh
python3 scripts/p337_direct_birth_arm_topology.py
python3 -m unittest tests.test_p337_direct_birth_arm_topology
```

The standard-library certificate writes `results/exact-direct-birth-arm-topology/latest.json`.  It reuses the existing exact subset homology oracle and does not run Monte Carlo.

Related issues: #337, #334.
