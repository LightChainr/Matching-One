# Digital Alexander rank oracle

Finite exact enumeration plus deterministic counterexample search; not a general-topology proof.

| geometry | N | configurations | joint `(r_b,r_w,q)` counts | weak failures | strong failures |
|---|---:|---:|---|---:|---:|
| axis-L2 | 4 | 16 | `r0_r2_q-1:7, r1_r1_q+0:4, r2_r0_q+1:5` | 0 | 0 |
| axis-L3 | 9 | 512 | `r0_r2_q-1:259, r1_r1_q+0:162, r2_r0_q+1:91` | 0 | 0 |
| gaussian-2-1 | 5 | 32 | `r0_r2_q-1:16, r1_r1_q+0:10, r2_r0_q+1:6` | 0 | 0 |
| diamond-L2 | 8 | 256 | `r0_r2_q-1:143, r1_r1_q+0:68, r2_r0_q+1:45` | 0 | 0 |
| c4-self-matching-3-1 | 10 | 1024 | `r0_r2_q-1:352, r1_r1_q+0:320, r2_r0_q+1:352` | 0 | 0 |

## Verdict

- common either/cross channel premise on every exhaustive configuration: `True`
- `2q = r_black-r_white` on every exhaustive configuration: `True`
- stronger `r_black+r_white=2` on every exhaustive configuration: `True`
- on every declared exhaustive geometry, the archived equality `q = r_black-r_white` occurs exactly when `q=0`, hence at rank pair `(1,1)`
- deterministic search: 85152 configurations on 160 integer-period tori; weak counterexamples: 0

The nine-case rank lemma proves that equality of the `either` and `cross` differences implies the weak rank identity. The finite oracle verifies the premise; it does not replace a digital Alexander-duality proof for arbitrary tori.

## Boundary

- The exhaustive checks establish finite configuration identities only on the declared quotients.
- The fixed-seed general-period scan is a counterexample search, not statistical evidence or a proof.
- A general result still requires a periodic digital Alexander/relative-homology argument with the 4/8 adjacency pair.
- No CFT field, universal amplitude, or closed form for square-site p_c is inferred.
