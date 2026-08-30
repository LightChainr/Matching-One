# Digital Alexander local bridge

Executable 16-face certificate plus the complementary-subsurface duality proof.

## Local 4/8 certificate

- face patterns checked: 16
- retained-diagonal black masks: [5, 10]
- redundant diagonals replaced inside one face: 6
- connectivity/replacement/embedding gate: `True`

A white diagonal is retained only when its endpoints are the only two white corners. In every other active case a white NN boundary path replaces it inside the same square, preserving its ambient homology class.

## Surface theorem

For complementary compact subsurfaces U,V of a closed oriented surface S, im(H1(V)->H1(S)) is the intersection-form orthogonal complement of im(H1(U)->H1(S)).

1. intersection pairing identifies the annihilator of im H1(U) with ker[H^1(S)->H^1(U)]
2. the long exact sequence of (S,U) identifies that kernel with the image of H^1(S,U)
3. excision gives H^1(S,U)=H^1(V,boundary V)
4. Poincare-Lefschetz duality gives H^1(V,boundary V)=H_1(V)
5. naturality identifies the resulting map with inclusion H_1(V)->H_1(S)

## Rank consequence

| r_black | r_white | q | weak residual |
|---:|---:|---:|---:|
| 0 | 2 | -1 | 0 |
| 1 | 1 | 0 | 0 |
| 2 | 0 | 1 | 0 |

Within the declared honest square-cell scope, `r_black+r_white=2` and therefore `2q=r_black-r_white`.

## Boundary

- The proof scope requires an honest periodic square-cell decomposition so the facewise regular-neighborhood construction is embedded.
- Tiny or short-period quotient degeneracies are covered by the separate finite oracle, not silently absorbed into the local-cell proof.
- The result identifies the lattice matching observable as an ambient-H1 rank imbalance; it does not identify a local CFT field or prove a selection rule for V_(2,2).
- No closed form for the square-site threshold follows from the homological balance identity.
